#include <Python.h>

//#define DEBUG

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "constants.h"

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>

#define GSL_EPS 0.5e-2
#define GSL_WSIZE 10000
#define GSL_KEY GSL_INTEG_GAUSS21

double gsl_epsilon=GSL_EPS;
int gsl_key=GSL_KEY;
int gsl_quiet=0;

double gsl_integ(double f(double),double x1,double x2,gsl_integration_workspace *w) {

  double res,err;
  gsl_function F;
  int r;

  // casts here are to stop the compiler complaining
  F.function=(double(*)(double, void *))f;
  F.params=NULL;
  //printf("gsl_integ called: %g %g\n",x1,x2);
  r=gsl_integration_qag(&F,x1,x2,0,gsl_epsilon,GSL_WSIZE,gsl_key,w,&res,&err);
  if (!gsl_quiet && r) fprintf(stderr,"gsl_integ error %i: returning %g, err is %g\n",r,res,err);
  return res;
}

double gsl_integ_i(double f(double),gsl_integration_workspace *w) {

  double res,err;
  gsl_function F;
  int r;

  F.function=(double(*)(double, void *))f;
  F.params=NULL;
  //printf("gsl_integ called: %g %g\n",x1,x2);
  r=gsl_integration_qagi(&F,0,GSL_EPS,GSL_WSIZE,w,&res,&err);
  if (!gsl_quiet && r) fprintf(stderr,"gsl_integ_i error %i: returning %g, err is %g\n",r,res,err);
  return res;
}

double gsl_integ_iu(double f(double),double x1,gsl_integration_workspace *w) {

  double res,err;
  gsl_function F;
  int r;

  F.function=(double(*)(double, void *))f;
  F.params=NULL;
  r=gsl_integration_qagiu(&F,x1,0,GSL_EPS,GSL_WSIZE,w,&res,&err);
  if (!gsl_quiet && r) fprintf(stderr,"gsl_integ_iu error %i: returning %g, err is %g\n",r,res,err);
  return res;
}

#define NUSE1 5
#define NUSE2 5
#define EPS 1.0e-10
#define FPMIN 1.0e-30
#define MAXIT 10000
#define XMIN 2.0

#define XF1 1.0e-4
#define XF2 0.22e+2
#define INFIN  1.0e2

#ifndef INTCO
#define INTCO (3*ECF)/(4*PI*M_EL*M_EL*M_EL*V_C*V_C*V_C*V_C)
#define EMCO (sqrt(3)*ECF*ECF*ECF)/(4.0*PI*EPS_0*V_C*M_EL)
#endif

int kp=0;
static double sinalpha;
static double axf1, axf2;
static double *fx;
static int xls=100;

double (*ne)(double);
double (*ng)(double);
double (*age_emax)(double);

double gmin,gmax,power,gbreak,gbvalue;
double EMIN,EMAX,EBREAK;
double n0_ext;
double BFIELD, nu;
double ebrd,gbrd;

double ageb,age;

gsl_integration_workspace *w1,*w2;

double ng_pow(double g) {

  if (g<gmin || g>gmax) return 0.0;
  else return pow(g,-power);
}

double ne_pow(double e) {

  /* ne(E) dE is the number of electrons with energy between E and E+dE */

  /* delta is the power-law index */
  /* normalise to 1 at some arbitrary energy */

  if ((e<EMIN) || (e>EMAX)) return 0.0;
  else return n0_ext*pow(e,-power);
}

double ne_age_jp(double e) {

  /* as above, but with a finite age single-burst model */
  /* loss is the term that corrects for loss. */
  /* JP version */

static double loss;

  if ((e<EMIN) || (e>EMAX)) return(0.0);
  else {
    loss=(4*THOMSON/(6*M_EL*M_EL*V_C*V_C*V_C*MU_0))*ageb*ageb*e*age;
    if (loss>=1.0) return(0.0);
    else return n0_ext*pow(e,-power)*pow((1.0-loss),power-2);
  }
}

double ne_age_kp(double e) {

  /* as above, but with a finite age single-burst model */
  /* loss is the term that corrects for loss. */
  /* KP version */

static double loss;

  if ((e<EMIN) || (e>EMAX)) return(0.0);
  else {
    loss=(THOMSON*sinalpha*sinalpha/(M_EL*M_EL*V_C*V_C*V_C*MU_0))*ageb*ageb*e*age;
    if (loss>=1.0) return(0.0);
    else return n0_ext*pow(e,-power)*pow((1.0-loss),power-2);
  }
}

double ng_age(double g) {

  /* as above, but with a finite age single-burst model */
  /* loss is the term that corrects for loss. */

static double loss;

  if (g<gmin || g>gmax) return 0.0;
  else {
    loss=(4*THOMSON/(6*M_EL*V_C*MU_0))*ageb*ageb*g*age;
    if (loss>=1.0) return(0.0);
    else return pow(g,-power)*pow((1.0-loss),power-2);
  }
}

double age_emax_jp(double limit) {

double em;

em=6.0*M_EL*M_EL*V_C*V_C*V_C*MU_0/(4*THOMSON*ageb*ageb*age);
if (em>limit) em=limit;
return(em);
}

double age_emax_kp(double limit) {

double em;

em=M_EL*M_EL*V_C*V_C*V_C*MU_0/(THOMSON*ageb*ageb*age*sinalpha*sinalpha);
if (em>limit) em=limit;
return(em);
}

double ng_break(double g) {
  if (g<gmin || g>gmax) return 0.0;
  else if (g>=gbreak) return gbrd*pow(g,-(power+gbvalue));
  else return pow(g,-power);
}

double ne_break(double e) {

  /* as above, but with a break */
  
  if ((e<EMIN) || (e>EMAX)) return 0.0;
  else if (e>=EBREAK) return n0_ext*ebrd*pow(e,-(power+gbvalue));
  else return n0_ext*pow(e,-power);
}

double ene(double e) {
  return e*ne(e);
}

double ffx(double x)

{
/* use the tabulated values of F(x) to return F(x) at any x. For x
   outside the table, use the asymptote of Pacholczyk (1970) for small
   x and 0 for large x. Otherwise return a log-linear two-point
   interpolation. */

  int i;
  double av,xv1,xv2;

  if (x>=XF2)
    return 0.0;
  else if (x<=XF1)
    return 4.0*PI*pow((0.5*x),(1.0/3.0))/(sqrt(3.0)*2.68357);
  /* this is a cheat because I'm having problems with the gamma fn */
  else {
    av=log(x);
    xv1=(av-axf1)*((double)xls-1.0)/(axf2-axf1);
    xv2=floor(xv1);
    i=xv2;
    xv1-=xv2;
    return exp(((1.0-xv1)*log(fx[i]))+(xv1*log(fx[i+1])));
  }
}

double ff(double x) {
/* given x, returns the modified Bessel function of order 5/3 */
  return gsl_sf_bessel_Knu(5.0/3.0,x);
}

void makefx(void) {

int i;
double xv, yv, av;

#ifdef DEBUG
printf("Making F(x) lookup table:\n"); 
#endif
axf1=log(XF1);
axf2=log(XF2);
for (i=0; i<xls; i++) {
  /* compute the x-value */
  av=axf1+((axf2-axf1)*((double)i)/(double)(xls-1));
  xv=exp(av);
  /* integrate ff from xv to infinity */
  yv=gsl_integ_iu(ff,xv,w1);
  //  if (yv==FAILTAG) break;
  /* once qromo1 has failed to converge, assume that the rest are 0 */
  /* F(x) = x times integral */
  fx[i]=xv*yv;
#ifdef DEBUG
  printf("%f %g\n",xv,fx[i]);
#endif
}

for (;i<xls; i++) fx[i]=0.0;
}

double integr(double e) {

/* the synchrotron emission integrand */
  return ne(e)*ffx(nu/(INTCO*BFIELD*sinalpha*e*e));
}

double synchpl(void) {

  static double emin,emax,res;

  emin=sqrt(nu/(BFIELD*sinalpha*INTCO*XF2));
  emax=3*sqrt(nu/(BFIELD*sinalpha*INTCO*XF1));
  if (emin<EMIN) emin=EMIN;
  if (emax>EMAX) emax=EMAX;
  if (age>0.0) emax=age_emax(emax);
  //  else emax=EMAX;
  /* printf("Using min energy %lg\n",emin); */
  /*  printf("in synchpl, emin = %g, emax = %g\n",emin,EMAX); */
  if (emin>emax) {
    //    printf("Synch. integration would be void; no suitable electrons (%g, %g).\n",emin,emax);
    res=0;
  } else
    res=gsl_integ(integr,emin,emax,w1);

/*  printf("returning %g\n",res); */
  return EMCO*BFIELD*sinalpha*res;
}

double synchin(double alpha) {

  double retval;

#ifdef DEBUG
  printf("Here i am in synchin, alpha is %g\n",alpha);
#endif
  sinalpha=sin(alpha);
  retval=sinalpha*synchpl();
#ifdef DEBUG
  printf("Retval is %g\n",retval);
#endif
  return retval;
}

double synchint(void) {

/* integrate over all pitch angles */

  return 0.5*gsl_integ(synchin,0.0,PI,w2);
}

double emiss_n(double n0, double b, double nu_arg) {

  /* similar to emiss_a, but the parameters of the electron
     distribution are fixed at init time, and can't be changed. */

  n0_ext=n0;
  BFIELD=b;
  nu=nu_arg;
#ifdef DEBUG
  printf("Calling synchint with %g %g %g\n",n0_ext,BFIELD,nu);
#endif
  return synchint();
}

/////////////////////// PYTHON API ///////////////////////////////

static PyObject *synch_setage(PyObject *self, PyObject *args) {
  if (!PyArg_ParseTuple(args, "dd", &age, &ageb))
    return NULL;
  return Py_BuildValue("d", age);
}

static PyObject *synch_setspectrum(PyObject *self, PyObject *args) {
  if (!PyArg_ParseTuple(args, "ddd", &gmin, &gmax, &power))
    return NULL;
  EMIN=gmin*M_EL*V_C*V_C;
  EMAX=gmax*M_EL*V_C*V_C;
  return Py_BuildValue("d", power);
}

static PyObject *synch_emiss(PyObject *self, PyObject *args) {
  const double norm, bfield, nu;
  double emiss;

  /* key function -- calculate emissivity */

  if (!PyArg_ParseTuple(args, "ddd", &norm, &bfield, &nu))
    return NULL;
#ifdef DEBUG
  printf("Calling emiss_n with values %g, %g, %g\n",norm,bfield,nu);
#endif
  emiss=emiss_n(norm,bfield,nu);
#ifdef DEBUG
  printf("Value returned was %g\n",emiss);
#endif
  return Py_BuildValue("d", emiss);
}

static PyObject *synch_ff(PyObject *self, PyObject *args) {

  double x;
  if (!PyArg_ParseTuple(args, "d", &x))
    return NULL;
  return Py_BuildValue("d", ff(x));
}

static PyObject *synch_fx(PyObject *self, PyObject *args) {

  double x;
  if (!PyArg_ParseTuple(args, "d", &x))
    return NULL;
  return Py_BuildValue("d", ffx(x));
}

static PyMethodDef SynchMethods[] = {
    {"ff",  synch_ff, METH_VARARGS,
     "Find modified Bessel fn of order 5/3."},
    {"fx",  synch_fx, METH_VARARGS,
     "Find F(x)."},
    {"emiss",  synch_emiss, METH_VARARGS,
     "Calculate an emissivity."},
    {"setage",  synch_setage, METH_VARARGS,
     "Set the synchrotron age to use (s) and ageing field (T)."},
    {"setspectrum",  synch_setspectrum, METH_VARARGS,
     "Set the gamma_min, gamma_max and power-law index."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC initsynch(void) {

  (void)Py_InitModule("synch", SynchMethods);

  // do some initialization

  w1=gsl_integration_workspace_alloc(GSL_WSIZE);
  w2=gsl_integration_workspace_alloc(GSL_WSIZE);
  gsl_set_error_handler_off(); // to avoid abort
  gsl_epsilon=1e-3;
  gsl_quiet=1;

  fx=calloc(xls,sizeof(double));
  #ifdef DEBUG
      printf("running makefx\n");
  #endif
  makefx();

  // these numeric parameters and the spectral function to be used should be changeable by the user

  gmin=10;
  gmax=100000;
  EMIN=gmin*M_EL*V_C*V_C;
  EMAX=gmax*M_EL*V_C*V_C;
  power=2.2;

  ne=ne_age_jp;
  age_emax=age_emax_jp;

  ageb=1e-9;
  age=1e7*365*86400;

}
