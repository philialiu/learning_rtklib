#ifndef EKF_H
#define EKF_H

#include <iostream>
#include <rtklib.h>

#define SQR(x)      ((x)*(x))
#define ERR_CBIAS   0.3         /* code bias error Std (m) */
#define MIN_EL      (5.0*D2R)   /* min elevation for measurement error (rad) */
#define SQRT(x)    ((x)<0.0||(x)!=(x)?0.0:sqrt(x))

/* --------- basic functions ----------- */
static void cmatd(double *matrx, int n);
static void cmatf(float *matrx, int n);
static void carrd(double *x, int n);
static void csat(obsd_t *obs, int n);
static int nextobsf(const obs_t *obs, int *i, int rcv);
char* charcat(const char* str1, const char* str2);
char* strpl(const char *str, const char *oldstr, const char *newstr);
static void ecef2lla(double *x, double *pos);

/* --------- spp --------- */
static double gettgd(int sat, const nav_t *nav, int type);
static double prange(const obsd_t *obs, const nav_t *nav, const prcopt_t *opt,
                     double *var);
static double varerr(const prcopt_t *opt, double el, int sys);
static int rescode(int iter, const obsd_t *obs, int n, const double *rs,
                   const double *dts, const double *vare, const int *svh,
                   const nav_t *nav, const double *x, const prcopt_t *opt,
                   double *v, double *H, double *var, double *azel, int *vsat,
                   double *resp, int *ns);

/* ------------ ekf ------------- */
static void ekfinit(sol_t *sol, ekfsol_t *esol);
static void predict(double tt, ekfsol_t *esol, prcopt_t *opt);
static int innocode(const obsd_t *obs, int n, const nav_t *nav, ekfsol_t *esol, 
                 const prcopt_t *opt, double *v, double *H, double *var);
static int filter_(const double *x, const double *P, const double *H,
                   const double *v, const double *R, int n, int m,
                   double *xp, double *Pp);
static void update(const obsd_t *obs, int n, const nav_t *nav,
                  const prcopt_t *opt, ekfsol_t *esol);

/*-------- file --------- */
static void outheader(FILE *fp, char **file, int n, const prcopt_t *popt,
                      const solopt_t *sopt, const obs_t *obss);
static int outhead(const char *outfile, char **infile, int n,
                   const prcopt_t *popt, const solopt_t *sopt, const obs_t *obs);
static FILE *openfile(const char *outfile);

#endif
