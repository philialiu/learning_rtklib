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
#include <Eigen/LU>

#define NX          (3*3+5*2)
#define SQR(x)      ((x)*(x))
#define ERR_CBIAS   0.3         /* code bias error Std (m) */
#define MIN_EL      (5.0*D2R)   /* min elevation for measurement error (rad) */

/* --------- basic functions ----------- */
static void cmatd(double *matrx, int n)
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

static void carrd(double *x, int n)
{
    for (int i=0;i<n;i++) {
        std::cout << x[i] << ' ';
    }
    std::cout << std::endl;
}

static void csat(obsd_t *obs, int n)
{
    int sat, sys, i;
    for (i=0;i<n;i++) {
        sat=obs[i].sat;
        sys=satsys(sat,NULL);
        std::cout << sys << std::endl;
    }
    
}

static int nextobsf(const obs_t *obs, int *i, int rcv)
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

char* charcat(const char* str1, const char* str2)
{
    char* str3 = new char[1024];
    strcpy(str3, str1);
    strcat(str3, str2);
    return str3;
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
    int i,j,nv=0,sat,sys,mask[NX-3]={0};
    
    trace(3,"resprng : n=%d\n",n);
    
    for (i=0;i<3;i++) rr[i]=x[i];
    dtr=x[3];
    
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
        
        if (iter>0) {
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
            H[j+nv*NX]=j<3?-e[j]:(j==3?1.0:0.0);
        }
        /* time system offset and receiver bias correction */
        if      (sys==SYS_GLO) {v[nv]-=x[4]; H[4+nv*NX]=1.0; mask[1]=1;}
        else if (sys==SYS_GAL) {v[nv]-=x[5]; H[5+nv*NX]=1.0; mask[2]=1;}
        else if (sys==SYS_CMP) {v[nv]-=x[6]; H[6+nv*NX]=1.0; mask[3]=1;}
        else if (sys==SYS_IRN) {v[nv]-=x[7]; H[7+nv*NX]=1.0; mask[4]=1;}
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
    /* constraint to avoid rank-deficient */
    for (i=0;i<NX-3;i++) {
        if (mask[i]) continue;
        v[nv]=0.0;
        for (j=0;j<NX;j++) H[j+nv*NX]=j==i+3?1.0:0.0;
        var[nv++]=0.01;
    }
    return nv;
}

/* ------------ ekf ------------- */
typedef struct {
    gtime_t time;       /* time (GPST) */
    double  x[NX];      /* position/velocity/acceleration/clock bias/clock drift */
    double  P[NX*NX];   /* covariance matrix */
} ekfsol_t;

/* spp result to ekf vars */
static void ekfinit(sol_t *sol, ekfsol_t *esol)
{
    int i, j;
    esol->time = sol->time;

    /* state */
    for (i=0;i<6;i++) {
        esol->x[i]=sol->rr[i]; // position/velocity
    }
    for (i=0;i<5;i++) {
        esol->x[i+9]=sol->dtr[i]; // clock bias/offset (s)
    }

    /* covariance matrix */
    for (i=0;i<NX*NX;i++) {
        esol->P[i] = 1000000;
    }
    for (i=0;i<3;i++) {
        esol->P[i+i*NX]=sol->qr[i];
    }
    esol->P[1]   =esol->P[NX]    =sol->qr[3];
    esol->P[NX+2]=esol->P[2*NX+1]=sol->qr[4];
    esol->P[2]   =esol->P[2*NX]  =sol->qr[5];

    // cmatd(esol->P, NX);
}

/* predict states */
static void predict(double tt, ekfsol_t *esol, prcopt_t *opt)
{
    double *F,*P,*FP,*x,*xp,x_[3],pos[3],Q[NX*NX],Qv[9];
    int i,j;

    /* state transition of position/velocity/clock bias/clock drift */
    F=eye(NX); P=mat(NX,NX); FP=mat(NX,NX); x=mat(NX,1); xp=mat(NX,1);
    
    for (i=0;i<6;i++) {
        F[i+(i+3)*NX]=tt;
    }
    for (i=0;i<3;i++) {
        F[i+(i+6)*NX]=SQR(tt)/2.0;
    }
    for (i=0;i<5;i++) {
        F[(i+3*3)+(i+3*3+5)*NX]=tt;
    }

    /* display F */
    // cmatd(F, NX);

    for (i=0;i<NX;i++) {
        x[i]=esol->x[i];
        for (j=0;j<NX;j++) {
            P[i+j*NX]=esol->P[i+j*NX];
        }
    }
    
    /* x=F*x, P=F*P*F+Q */
    matmul("NN",NX,1,NX,1.0,F,x,0.0,xp);
    matmul("NN",NX,NX,NX,1.0,F,P,0.0,FP);
    matmul("NT",NX,NX,NX,1.0,FP,F,0.0,P);
    
    for (i=0;i<NX;i++) {
        esol->x[i]=xp[i];
        for (j=0;j<NX;j++) {
            esol->P[i+j*NX]=P[i+j*NX];
        }
    }

    // /* process noise added to only acceleration */
    Q[0]=Q[4]=SQR(opt->prn[3])*fabs(tt);
    Q[8]=SQR(opt->prn[4])*fabs(tt);
    for (i=0;i<3;i++) {
        x_[i] = xp[i];
    }
    ecef2pos(x_,pos);
    covecef(pos,Q,Qv);
    for (i=0;i<3;i++) for (j=0;j<3;j++) {
        esol->P[i+6+(j+6)*NX]+=Qv[i+j*3];
    }

    free(F); free(P); free(FP); free(x); free(xp);
    carrd(esol->x, NX);
}

/* pseudorange residual */
static void inno(const obsd_t *obs, int n, const nav_t *nav, ekfsol_t *esol, 
                 const prcopt_t *opt, double *v, double *H)
{
    int i,nv,ns,svh[MAXOBS],vsat[MAXOBS]={0};
    double *rs,*dts,*var,*azel,*resp;
    double *v,*H;
    double x[3+5]={0};

    for (i=0;i<3;i++) {
        x[i] = esol->x[i];
    }
    for (i=0;i<5;i++) {
        x[i+3] = esol->x[i+9];
    }

    satposs(esol->time,obs,n,nav,opt->sateph,rs,dts,var,svh);

    azel=zeros(2,n); v=mat(n+4,1); H=mat(NX,n+4);
    nv=rescode(i,obs,n,rs,dts,var,svh,nav,x,opt,v,H,var,azel,vsat,resp,
                   &ns);

    free(rs); free(dts); free(var); free(azel); free(resp);
}

// TODO
static int update(const obsd_t *obs, int n, const double *rs, const double *dts,
                  const double *vare, const int *svh, const nav_t *nav,
                  const prcopt_t *opt, ekfsol_t *esol, double *azel, int *vsat,
                  double *resp, char *msg)
{
    double x[NX]={0},dx[NX],Q[NX*NX],*v,*H,*var,sig;
    v=mat(n+4,1); H=mat(NX,n+4); var=mat(n+4,1);

    inno(obs, n, nav, esol, opt, v, H);
}

int main(int argc, char **argv)
{
    gtime_t t0 = {0}, ts = {0}, te = {0};

    char filepath[] = "/home/philia/Documents/nav_data/uav_20211211/gnc01/data_20211211_150700/processed/";
    char file1[] = "eph_202112110655.nav";
    char file2[] = "ublox_20211211_14.obs";

    obs_t obs = {0};
    nav_t nav = {0};
    sta_t sta = {""};

    prcopt_t opt = prcopt_default;
    opt.navsys = SYS_GPS|SYS_GLO|SYS_GAL|SYS_CMP|SYS_QZS|SYS_IRN; // use all satellite systems
    sol_t sol;
    ekfsol_t esol;
    char msg[128];

    int m = 0;

    int stat1 = readrnxt(charcat(filepath, file1), 1, ts, te, 0.0, "", &obs, &nav, &sta);
    int stat2 = readrnxt(charcat(filepath, file2), 1, t0, t0, 0.0, "", &obs, &nav, &sta);

    if (stat1 && stat2)
    {
        std::cout << "Files are successfully opened!" << std::endl;
    }
    else
    {
        std::cout << "Fail to open files." << std::endl;
    }

    /* first epoch: spp */
    int i = 0;
    int n = nextobsf(&obs, &i, 1);
    int ret = pntpos(&obs.data[i], n, &nav, &opt, &sol, NULL, NULL, msg);
    double tt = timediff(obs.data[i + n].time, obs.data[i].time);
    ekfinit(&sol, &esol);
    if (ret == 1) // 1：OK, 0: error
    {
        double ep[6] = {0};
        time2epoch(sol.time, ep);
        printf("%.0lf,%.0lf,%.0lf,%.0lf,%.0lf,%.0lf,%lf,%lf,%lf,%lf,%lf,%lf,\n", ep[0], ep[1], ep[2], ep[3], ep[4], ep[5],
            sol.rr[0], sol.rr[1], sol.rr[2], sol.rr[3], sol.rr[4], sol.rr[5]);
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
        // update();
    }

    return 0;
}