/*------------------------------------------------------------------------------
* ekf.c : spp ekf
*
*          Copyright (C) 2022 by Philia Liu, All rights reserved.
*
* reference :
*     [1] 
*
*-----------------------------------------------------------------------------*/
#include <iostream>
#include <rtklib.h>
#include <fstream>

#define NX          (3*3+5)
#define SQR(x)      ((x)*(x))
#define ERR_CBIAS   0.3         /* code bias error Std (m) */
#define MIN_EL      (5.0*D2R)   /* min elevation for measurement error (rad) */

typedef struct {
    sol_t   sol;        /* solution */
    double  x[14];      /* position/velocity/acceleration/clock bias/clock drift */
    double  P[14*14];   /* covariance matrix */
} ekfsol_t;

/* --------- basic functions ----------- */
static void cmatd(double *matrx, int n) // output matrix
{
    for (int i=0;i<n;i++) {
        for (int j=0;j<n;j++) {
            std::cout << matrx[i*n+j] << ' ';
        }
        std::cout << std::endl;
    }
}

static void cmatf(float *matrx, int n)
{
    for (int i=0;i<n;i++) {
        for (int j=0;j<n;j++) {
            std::cout << matrx[i*n+j] << ' ';
        }
        std::cout << std::endl;
    }
}

static void carrd(double *x, int n) // output array
{
    for (int i=0;i<n;i++) {
        std::cout << x[i] << ' ';
    }
    std::cout << std::endl;
}

static void csat(obsd_t *obs, int n) // output satellite system
{
    int sat, sys, i;
    for (i=0;i<n;i++) {
        sat=obs[i].sat;
        sys=satsys(sat,NULL);
        std::cout << sys << std::endl;
    }
}

static int nextobsf(const obs_t *obs, int *i, int rcv) // search next observation data index
{
    double tt;
    int n;

    for (; *i < obs->n; (*i)++)
        if (obs->data[*i].rcv == rcv)
            break;
    for (n = 0; *i + n < obs->n; n++)
    {
        /* time difference (t1-t2) (s) */
        tt = timediff(obs->data[*i + n].time, obs->data[*i].time);
        if (obs->data[*i + n].rcv != rcv || tt > DTTOL)
            break;
    }
    return n;
}

char* charcat(const char* str1, const char* str2) // concatenate string
{
    char* str3 = new char[1024];
    strcpy(str3, str1);
    strcat(str3, str2);
    return str3;
}

static void ecef2lla(double *x, double *pos)
{
    double xx[3];
    for (int j=0;j<3;j++) xx[j]=x[j];
    ecef2pos(x,pos);
    pos[0] *= R2D;
    pos[1] *= R2D;
}

/* --------- spp --------- */
static double gettgd(int sat, const nav_t *nav, int type)
{
    int i,sys=satsys(sat,NULL);
    
    if (sys==SYS_GLO) {
        for (i=0;i<nav->ng;i++) {
            if (nav->geph[i].sat==sat) break;
        }
        return (i>=nav->ng)?0.0:-nav->geph[i].dtaun*CLIGHT;
    }
    else {
        for (i=0;i<nav->n;i++) {
            if (nav->eph[i].sat==sat) break;
        }
        return (i>=nav->n)?0.0:nav->eph[i].tgd[type]*CLIGHT;
    }
}

static double prange(const obsd_t *obs, const nav_t *nav, const prcopt_t *opt,
                     double *var)
{
    double P1,P2,gamma,b1,b2;
    int sat,sys;
    
    sat=obs->sat;
    sys=satsys(sat,NULL);
    P1=obs->P[0];
    P2=obs->P[1];
    *var=0.0;
    
    if (P1==0.0||(opt->ionoopt==IONOOPT_IFLC&&P2==0.0)) return 0.0;
    
    /* P1-C1,P2-C2 DCB correction */
    if (sys==SYS_GPS||sys==SYS_GLO) {
        if (obs->code[0]==CODE_L1C) P1+=nav->cbias[sat-1][1]; /* C1->P1 */
        if (obs->code[1]==CODE_L2C) P2+=nav->cbias[sat-1][2]; /* C2->P2 */
    }
    if (opt->ionoopt==IONOOPT_IFLC) { /* dual-frequency */
        
        if (sys==SYS_GPS||sys==SYS_QZS) { /* L1-L2,G1-G2 */
            gamma=SQR(FREQ1/FREQ2);
            return (P2-gamma*P1)/(1.0-gamma);
        }
        else if (sys==SYS_GLO) { /* G1-G2 */
            gamma=SQR(FREQ1_GLO/FREQ2_GLO);
            return (P2-gamma*P1)/(1.0-gamma);
        }
        else if (sys==SYS_GAL) { /* E1-E5b */
            gamma=SQR(FREQ1/FREQ7);
            if (getseleph(SYS_GAL)) { /* F/NAV */
                P2-=gettgd(sat,nav,0)-gettgd(sat,nav,1); /* BGD_E5aE5b */
            }
            return (P2-gamma*P1)/(1.0-gamma);
        }
        else if (sys==SYS_CMP) { /* B1-B2 */
            gamma=SQR(((obs->code[0]==CODE_L2I)?FREQ1_CMP:FREQ1)/FREQ2_CMP);
            if      (obs->code[0]==CODE_L2I) b1=gettgd(sat,nav,0); /* TGD_B1I */
            else if (obs->code[0]==CODE_L1P) b1=gettgd(sat,nav,2); /* TGD_B1Cp */
            else b1=gettgd(sat,nav,2)+gettgd(sat,nav,4); /* TGD_B1Cp+ISC_B1Cd */
            b2=gettgd(sat,nav,1); /* TGD_B2I/B2bI (m) */
            return ((P2-gamma*P1)-(b2-gamma*b1))/(1.0-gamma);
        }
        else if (sys==SYS_IRN) { /* L5-S */
            gamma=SQR(FREQ5/FREQ9);
            return (P2-gamma*P1)/(1.0-gamma);
        }
    }
    else { /* single-freq (L1/E1/B1) */
        *var=SQR(ERR_CBIAS);
        
        if (sys==SYS_GPS||sys==SYS_QZS) { /* L1 */
            b1=gettgd(sat,nav,0); /* TGD (m) */
            return P1-b1;
        }
        else if (sys==SYS_GLO) { /* G1 */
            gamma=SQR(FREQ1_GLO/FREQ2_GLO);
            b1=gettgd(sat,nav,0); /* -dtaun (m) */
            return P1-b1/(gamma-1.0);
        }
        else if (sys==SYS_GAL) { /* E1 */
            if (getseleph(SYS_GAL)) b1=gettgd(sat,nav,0); /* BGD_E1E5a */
            else                    b1=gettgd(sat,nav,1); /* BGD_E1E5b */
            return P1-b1;
        }
        else if (sys==SYS_CMP) { /* B1I/B1Cp/B1Cd */
            if      (obs->code[0]==CODE_L2I) b1=gettgd(sat,nav,0); /* TGD_B1I */
            else if (obs->code[0]==CODE_L1P) b1=gettgd(sat,nav,2); /* TGD_B1Cp */
            else b1=gettgd(sat,nav,2)+gettgd(sat,nav,4); /* TGD_B1Cp+ISC_B1Cd */
            return P1-b1;
        }
        else if (sys==SYS_IRN) { /* L5 */
            gamma=SQR(FREQ9/FREQ5);
            b1=gettgd(sat,nav,0); /* TGD (m) */
            return P1-gamma*b1;
        }
    }
    return P1;
}

static double varerr(const prcopt_t *opt, double el, int sys)
{
    double fact,varr;
    fact=sys==SYS_GLO?EFACT_GLO:(sys==SYS_SBS?EFACT_SBS:EFACT_GPS);
    if (el<MIN_EL) el=MIN_EL;
    varr=SQR(opt->err[0])*(SQR(opt->err[1])+SQR(opt->err[2])/sin(el));
    if (opt->ionoopt==IONOOPT_IFLC) varr*=SQR(3.0); /* iono-free */
    return SQR(fact)*varr;
}

static int rescode(int iter, const obsd_t *obs, int n, const double *rs,
                   const double *dts, const double *vare, const int *svh,
                   const nav_t *nav, const double *x, const prcopt_t *opt,
                   double *v, double *H, double *var, double *azel, int *vsat,
                   double *resp, int *ns)
{
    gtime_t time;
    double r,freq,dion,dtrp,vmeas,vion,vtrp,rr[3],pos[3],dtr,e[3],P;
    int i,j,nv=0,sat,sys,mask[NX-9]={0};
    
    trace(3,"resprng : n=%d\n",n);
    
    for (i=0;i<3;i++) rr[i]=x[i];
    dtr=x[9];
    
    ecef2pos(rr,pos);
    
    for (i=*ns=0;i<n&&i<MAXOBS;i++) {
        vsat[i]=0; azel[i*2]=azel[1+i*2]=resp[i]=0.0;
        time=obs[i].time;
        sat=obs[i].sat;
        if (!(sys=satsys(sat,NULL))) continue;
        
        /* reject duplicated observation data */
        if (i<n-1&&i<MAXOBS-1&&sat==obs[i+1].sat) {
            trace(2,"duplicated obs data %s sat=%d\n",time_str(time,3),sat);
            i++;
            continue;
        }
        /* excluded satellite? */
        if (satexclude(sat,vare[i],svh[i],opt)) continue;
        
        /* geometric distance */
        if ((r=geodist(rs+i*6,rr,e))<=0.0) continue;
        
        if (iter>=0) {
            /* test elevation mask */
            if (satazel(pos,e,azel+i*2)<opt->elmin) continue;
            
            /* ionospheric correction */
            if (!ionocorr(time,nav,sat,pos,azel+i*2,opt->ionoopt,&dion,&vion)) {
                continue;
            }
            if ((freq=sat2freq(sat,obs[i].code[0],nav))==0.0) continue;
            dion*=SQR(FREQ1/freq);
            vion*=SQR(FREQ1/freq);
            
            /* tropospheric correction */
            if (!tropcorr(time,nav,pos,azel+i*2,opt->tropopt,&dtrp,&vtrp)) {
                continue;
            }
        }
        /* psendorange with code bias correction */
        if ((P=prange(obs+i,nav,opt,&vmeas))==0.0) continue;
        
        /* pseudorange residual */
        v[nv]=P-(r+dtr-CLIGHT*dts[i*2]+dion+dtrp);

        /* design matrix */
        for (j=0;j<NX;j++) {
            H[j+nv*NX]=j<3?-e[j]:(j==9?1.0:0.0);
        }

        /* time system offset and receiver bias correction */
        if      (sys==SYS_GLO) {v[nv]-=x[10]; H[10+nv*NX]=1.0; mask[1]=1;}
        else if (sys==SYS_GAL) {v[nv]-=x[11]; H[11+nv*NX]=1.0; mask[2]=1;}
        else if (sys==SYS_CMP) {v[nv]-=x[12]; H[12+nv*NX]=1.0; mask[3]=1;}
        else if (sys==SYS_IRN) {v[nv]-=x[13]; H[13+nv*NX]=1.0; mask[4]=1;}
#if 0 /* enable QZS-GPS time offset estimation */
        else if (sys==SYS_QZS) {v[nv]-=x[8]; H[8+nv*NX]=1.0; mask[5]=1;}
#endif
        else mask[0]=1;
        
        vsat[i]=1; resp[i]=v[nv]; (*ns)++;
        
        /* variance of pseudorange error */
        var[nv++]=varerr(opt,azel[1+i*2],sys)+vare[i]+vmeas+vion+vtrp;
        
        trace(4,"sat=%2d azel=%5.1f %4.1f res=%7.3f sig=%5.3f\n",obs[i].sat,
              azel[i*2]*R2D,azel[1+i*2]*R2D,resp[i],sqrt(var[nv-1]));
    }

    // test
    // for (int k=0;k<(n+4);k++) {
    //     for (j=0;j<NX;j++) {
    //         printf("%lf ", H[k*NX+j]); 
    //     }
    //     printf("\n");
    // }
    // printf("----\n");
    // for (int k=0;k<nv;k++) {
    //     printf("%lf ", v[k]);
    // }
    // printf("\n----\n");

    /* constraint to avoid rank-deficient */
    for (i=0;i<NX-9;i++) {
        if (mask[i]) continue;
        v[nv]=0.0;
        for (j=0;j<NX;j++) H[j+nv*NX]=j==i+9?1.0:0.0;
        var[nv++]=0.01;
    }

    // test
    // for (int k=0;k<(n+4);k++) {
    //     for (j=0;j<NX;j++) {
    //         printf("%lf ", H[k*NX+j]); 
    //     }
    //     printf("\n");
    // }
    // printf("----\n");
    // for (int k=0;k<nv;k++) {
    //     printf("%lf ", v[k]);
    // }
    // printf("\n----\n");

    return nv;
}

/* ------------ ekf ------------- */
/* spp result to ekf vars */
static void ekfinit(sol_t *sol, ekfsol_t *esol)
{
    int i,j;
    double *P;
    esol->sol.time = sol->time;

    /* state */
    for (i=0;i<6;i++) {
        esol->x[i]=sol->rr[i]; // position/velocity
    }
    for (i=0;i<5;i++) {
        esol->x[i+9]=sol->dtr[i]*CLIGHT; // clock bias/offset (m)
    }

    /* covariance matrix */
    P=zeros(NX,NX);
    for (i=0;i<3;i++) {
        P[i+i*NX]=sol->qr[i];
    }
    P[1]   =P[NX]    =sol->qr[3];
    P[NX+2]=P[2*NX+1]=sol->qr[4];
    P[2]   =P[2*NX]  =sol->qr[5];
    for (i=3;i<NX;i++) {
        P[i+i*NX] = 1;
    }
    matcpy(esol->P,P,NX,NX);
    free(P);

    // cmatd(esol->P, NX);
}

/* predict states */
static void predict(double tt, ekfsol_t *esol, prcopt_t *opt)
{
    double *F,*P,*FP,*x,*xp,x_[3],pos[3],Q[9]={0},Qv[9];
    int i,j;

    /* state transition of position/velocity/clock bias/clock drift */
    F=eye(NX); P=mat(NX,NX); FP=mat(NX,NX); x=mat(NX,1); xp=mat(NX,1);
    
    for (i=0;i<6;i++) {
        F[i+(i+3)*NX]=tt;
    }
    for (i=0;i<3;i++) {
        F[i+(i+6)*NX]=SQR(tt)/2.0;
    }

    /* display F */
    // cmatd(F, NX);

    matcpy(P, esol->P, NX, NX);
    // cmatd(P, NX);
    
    /* x=F*x, P=F*P*F+Q */
    matmul("NN",NX,1,NX,1.0,F,x,0.0,xp);
    matmul("NN",NX,NX,NX,1.0,F,P,0.0,FP);
    matmul("NT",NX,NX,NX,1.0,FP,F,0.0,P);

    // cmatd(P, NX);

    /* process noise added to only acceleration */
    Q[0]=Q[4]=SQR(opt->prn[3])*fabs(tt);
    Q[8]=SQR(opt->prn[4])*fabs(tt);
    for (i=0;i<3;i++) {
        x_[i] = xp[i];
    }
    ecef2pos(x_,pos);
    covecef(pos,Q,Qv);
    for (i=0;i<3;i++) for (j=0;j<3;j++) {
        P[i+6+(j+6)*NX]+=Qv[i+j*3];
    }

    matcpy(esol->P, P, NX, NX);

    // test
    // cmatd(esol->P, NX);

    free(F); free(P); free(FP); free(x); free(xp);
    // carrd(esol->x, NX);
}

/* pseudorange residual */
static int inno(const obsd_t *obs, int n, const nav_t *nav, ekfsol_t *esol, 
                 const prcopt_t *opt, double *v, double *H, double *var)
{
    esol->sol.time=obs[0].time;
    int i,nv,ns,svh[MAXOBS],vsat[MAXOBS]={0};
    double *rs,*dts,*azel,*resp,*vare;
    double x[NX]={0};

    for (i=0;i<NX;i++) {
        x[i] = esol->x[i];
    }

    rs=mat(6,n); dts=mat(2,n); vare=mat(1,n);
    satposs(esol->sol.time,obs,n,nav,opt->sateph,rs,dts,vare,svh);

    azel=zeros(2,n); resp=mat(1,n);
    nv=rescode(0,obs,n,rs,dts,vare,svh,nav,x,opt,v,H,var,azel,vsat,resp,
                   &ns);

    free(rs); free(dts); free(azel); free(resp);
    return nv;
}

static int filter_(const double *x, const double *P, const double *H,
                   const double *v, const double *R, int n, int m,
                   double *xp, double *Pp)
{
    double *F=mat(n,m),*Q=mat(m,m),*K=mat(n,m),*I=eye(n);
    int info;
    
    matcpy(Q,R,m,m);
    matcpy(xp,x,n,1);
    matmul("NN",n,m,n,1.0,P,H,0.0,F);       /* Q=H'*P*H+R */
    matmul("TN",m,m,n,1.0,H,F,1.0,Q);

    // test
    // for (int i=0;i<n;i++) {
    //     for (int j=0;j<n;j++) {
    //         printf("%lf ", P[i*n+j]); 
    //     }
    //     printf("\n");
    // }
    // printf("----\n");

    // for (int i=0;i<m;i++) {
    //     for (int j=0;j<n;j++) {
    //         printf("%lf ", H[i*n+j]); 
    //     }
    //     printf("\n");
    // }
    // printf("----\n");

    // for (int i=0;i<m;i++) {
    //     for (int j=0;j<m;j++) {
    //         printf("%lf ", Q[i*m+j]); 
    //     }
    //     printf("\n");
    // }
    // printf("----\n");

    if (!(info=matinv(Q,m))) {
        matmul("NN",n,m,m,1.0,F,Q,0.0,K);   /* K=P*H*Q^-1 */
        matmul("NN",n,1,m,1.0,K,v,1.0,xp);  /* xp=x+K*v */
        matmul("NT",n,n,m,-1.0,K,H,1.0,I);  /* Pp=(I-K*H')*P */
        matmul("NN",n,n,n,1.0,I,P,0.0,Pp);
    }
    free(F); free(Q); free(K); free(I);
    return info;
}

/* measurement update */
static void update(const obsd_t *obs, int n, const nav_t *nav,
                  const prcopt_t *opt, ekfsol_t *esol)
{
    gtime_t time=obs[0].time;
    int info,i,j;
    double x[NX]={0},P[NX*NX]={0},*v,*H,*R,*var,*v_,*H_,x_[NX]={0},P_[NX*NX]={0};
    v=zeros(n+4,1); H=zeros(NX,n+4); var=zeros(n+4,1);

    /* innovation */
    int nv=inno(obs, n, nav, esol, opt, v, H, var);

    /* measurement noise */
    R=zeros(nv,nv);
    for (i=0;i<nv;i++) {
        R[i+i*nv]=var[i];
    }

    // cmatd(R, nv);

    /* prune matrix */
    v_=zeros(nv,1); H_=zeros(NX,nv);
    for (i=0;i<nv;i++) {
        v_[i]=v[i];
    }
    for (i=0;i<NX*nv;i++) {
        H_[i]=H[i];
    }
    free(v); free(H);

    /* kalman filter */
    matcpy(x,esol->x,NX,1);
    matcpy(P,esol->P,NX,NX);
    // cmatd(P, NX);
    // std::cout << "....." << std::endl;
    if ((info=filter_(x,P,H_,v_,R,NX,nv,x_,P_))) {
        std::cout << "kf update error" << std::endl;
    }
    // cmatd(P, NX);
    std::cout << time.time << ", P=" << P[0] << std::endl;

    /* update state and covariance matrix */
    matcpy(esol->x,x_,NX,1);
    matcpy(esol->P,P_,NX,NX);
    free(v_); free(H_); free(R); free(var);
    esol->sol.time=time;
    for (i=0;i<6;i++) esol->sol.rr[i]=esol->x[i];
    for (i=0;i<3;i++) esol->sol.qr[i]=esol->P[i+i*NX];
    esol->sol.qr[3]=esol->P[1];
    esol->sol.qr[4]=esol->P[NX+2];
    esol->sol.qr[5]=esol->P[2];
    for (i=0;i<3;i++) esol->sol.qv[i]=esol->P[(i+3)+(i+3)*NX];
    esol->sol.qv[3]=esol->P[NX*3+4];
    esol->sol.qv[4]=esol->P[NX*4+5];
    esol->sol.qv[5]=esol->P[NX*3+5];
    for (i=0;i<5;i++) esol->sol.dtr[i]=esol->x[i+9];
}

/*-------- file --------- */
static void outheader(FILE *fp, char **file, int n, const prcopt_t *popt,
                      const solopt_t *sopt, const obs_t *obss)
{
    const char *s1[]={"GPST","UTC","JST"};
    gtime_t ts,te;
    double t1,t2;
    int i,j,w1,w2;
    char s2[32],s3[32];
    
    trace(3,"outheader: n=%d\n",n);
    
    if (sopt->posf==SOLF_NMEA||sopt->posf==SOLF_STAT) {
        return;
    }
    if (sopt->outhead) {
        if (!*sopt->prog) {
            fprintf(fp,"%s program   : RTKLIB ver.%s\n",COMMENTH,VER_RTKLIB);
        }
        else {
            fprintf(fp,"%s program   : %s\n",COMMENTH,sopt->prog);
        }
        for (i=0;i<n;i++) {
            fprintf(fp,"%s inp file  : %s\n",COMMENTH,file[i]);
        }
        for (i=0;i<obss->n;i++)    if (obss->data[i].rcv==1) break;
        for (j=obss->n-1;j>=0;j--) if (obss->data[j].rcv==1) break;
        if (j<i) {fprintf(fp,"\n%s no rover obs data\n",COMMENTH); return;}
        ts=obss->data[i].time;
        te=obss->data[j].time;
        t1=time2gpst(ts,&w1);
        t2=time2gpst(te,&w2);
        if (sopt->times>=1) ts=gpst2utc(ts);
        if (sopt->times>=1) te=gpst2utc(te);
        if (sopt->times==2) ts=timeadd(ts,9*3600.0);
        if (sopt->times==2) te=timeadd(te,9*3600.0);
        time2str(ts,s2,1);
        time2str(te,s3,1);
        fprintf(fp,"%s obs start : %s %s (week%04d %8.1fs)\n",COMMENTH,s2,s1[sopt->times],w1,t1);
        fprintf(fp,"%s obs end   : %s %s (week%04d %8.1fs)\n",COMMENTH,s3,s1[sopt->times],w2,t2);
    }
    outprcopt(fp,popt);
    if (sopt->outhead||sopt->outopt) fprintf(fp,"%s\n",COMMENTH);
    
    outsolhead(fp,sopt);
}

static int outhead(const char *outfile, char **infile, int n,
                   const prcopt_t *popt, const solopt_t *sopt, const obs_t *obs)
{
    FILE *fp=stdout;
    
    trace(3,"outhead: outfile=%s n=%d\n",outfile,n);
    
    if (*outfile) {
        createdir(outfile);
        
        if (!(fp=fopen(outfile,"wb"))) {
            std::cout << "error : cannot open output file" << std::endl;
            return 0;
        }
    }
    /* output header */
    outheader(fp,infile,n,popt,sopt,obs);
    
    if (*outfile) fclose(fp);
    
    return 1;
}

static FILE *openfile(const char *outfile)
{
    trace(3,"openfile: outfile=%s\n",outfile);
    
    return !*outfile?stdout:fopen(outfile,"ab");
}


int main(int argc, char **argv)
{
    gtime_t t0 = {0}, ts = {0}, te = {0};

    char filepath[] = "/home/philia/Documents/nav_data/uav_20211211/gnc01/data_20211211_150700/processed/";
    char filename1[] = "eph_202112110655.nav";
    char filename2[] = "ublox_20211211_14.obs";
    char outfilename[] = "ublox_20211211_14_kf.pos";
    char *file1_ = charcat(filepath, filename1);
    char *file2_ = charcat(filepath, filename2);
    char *infile[2];
    char *outfile=charcat(filepath, outfilename);
    infile[0] = &file1_[0]; infile[1] = &file2_[0];
    FILE *fp;
    fp=openfile(outfile);

    obs_t obs = {0};
    nav_t nav = {0};
    sta_t sta = {""};

    prcopt_t opt = prcopt_default; // positioning option
    solopt_t solopt = solopt_default; // solution output options
    opt.navsys = SYS_GPS|SYS_GLO|SYS_GAL|SYS_CMP|SYS_QZS|SYS_IRN; // use all satellite systems
    opt.ionoopt = 1;
    opt.tropopt = 1;
    
    sol_t sol;
    ekfsol_t esol;
    char msg[128];

    int m = 0;
    double tt, rb[3]={0};

    int stat1 = readrnxt(file1_, 1, ts, te, 0.0, "", &obs, &nav, &sta);
    int stat2 = readrnxt(file2_, 1, t0, t0, 0.0, "", &obs, &nav, &sta);

    if (stat1 && stat2)
    {
        std::cout << "Files are successfully opened!" << std::endl;
    }
    else
    {
        std::cout << "Fail to open files." << std::endl;
    }

    /* write header to output file */
    outhead(outfile,infile,2,&opt,&solopt,&obs);

    /* first epoch: spp */
    int i = 0;
    int n = nextobsf(&obs, &i, 1);
    int ret = pntpos(&obs.data[i], n, &nav, &opt, &sol, NULL, NULL, msg);
    outsol(fp,&sol,rb,&solopt);
    tt = timediff(obs.data[i + n].time, obs.data[i].time);
    ekfinit(&sol, &esol);
    if (ret == 1) // 1：OK, 0: error
    {
        double ep[6] = {0},pos[3];
        time2epoch(sol.time, ep);
        ecef2lla(sol.rr,pos);
        printf("%.0lf,%.0lf,%.0lf,%.0lf,%.0lf,%.0lf,%lf,%lf,%lf,%lf,%lf,%lf,\n", ep[0], ep[1], ep[2], ep[3], ep[4], ep[5],
            pos[0], pos[1], pos[2], sol.rr[3], sol.rr[4], sol.rr[5]);
    }
    else
    {
        printf("ret: %d, msg:%s\n", ret, msg);
    }

    // test
    // csat(&obs.data[i], n);

    /* ekf */
    for (int i = n; (m = nextobsf(&obs, &i, 1)) > 0; i += m)
    {
        predict(tt, &esol, &opt);
        update(&obs.data[i], m, &nav, &opt, &esol);
        outsol(fp, &esol.sol, rb, &solopt);
        tt = timediff(obs.data[i + m].time, obs.data[i].time);
        double ep[6]={0},pos[3];
        ecef2lla(esol.x,pos);
        time2epoch(esol.sol.time, ep);
        printf("%.0lf,%.0lf,%.0lf,%.0lf,%.0lf,%.0lf,%lf,%lf,%lf,%lf,%lf,%lf,\n", ep[0], ep[1], ep[2], ep[3], ep[4], ep[5],
            pos[0], pos[1], pos[2], esol.x[3], esol.x[4], esol.x[5]);
    }
    fclose(fp);

    return 0;
}