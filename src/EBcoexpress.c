#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

SEXP f0worker(SEXP zvalsx, SEXP sizes3x, SEXP thetax, SEXP dx)
{
  int d = asInteger(dx);
  double *zvals = REAL(zvalsx);
  double *ivs = REAL(sizes3x);
  double *theta = REAL(thetax);
  double mu = theta[0];
  double tau2 = theta[1]*theta[1];
  double tot = 0;
  int i;

  double bit1 = 0;
  for(i=0; i < d; i++)
    bit1 = bit1 + ivs[i];
  bit1 = 1 + tau2*bit1;
  
  double prodvs = 1;
  for(i=0; i < d; i++)
    prodvs = prodvs/ivs[i];

  double bit2a = 0;
  for(i=0; i < d; i++)
    bit2a = bit2a + (zvals[i]-mu)*(zvals[i]-mu)*ivs[i];

  double bit2b = 0;
  for(i=0; i < d; i++)
    bit2b = bit2b + (zvals[i]-mu)*ivs[i];
  bit2b = (bit2b*bit2b*tau2)/bit1;

  tot = tot - d*log(M_2PI)/2;
  tot = tot - log(prodvs*bit1)/2;
  tot = tot - (bit2a-bit2b)/2;

  return ScalarReal(tot);
}

SEXP bwmcCworker(SEXP Xx, SEXP rnx, SEXP cnx, SEXP medsx, SEXP madsx)
{
  double *X = REAL(Xx);
  int rn = asInteger(rnx);
  int cn = asInteger(cnx);
  double *meds = REAL(medsx);
  double *mads = REAL(madsx);

  int i,j;
  SEXP XMs, Us, As, bidiag, bicor;
  double *rXMs, *rUs, *rAs, *rbidiag, *rbicor, temp, top, bot, tot1, tot2, s1, s2;
  PROTECT(XMs = allocMatrix(REALSXP,rn,cn));
  PROTECT(Us = allocMatrix(REALSXP,rn,cn));
  PROTECT(As = allocMatrix(REALSXP,rn,cn));
  PROTECT(bidiag = allocVector(REALSXP,cn));
  PROTECT(bicor = allocMatrix(REALSXP,cn,cn));
  
  rXMs = REAL(XMs);
  rUs = REAL(Us);
  rAs = REAL(As);
  rbidiag = REAL(bidiag);
  rbicor = REAL(bicor);
  
  for (j=0; j < cn; j++) {
    for(i=0; i < rn; i++) {
      temp = X[i+ rn*j] - meds[j];
      rXMs[i + rn*j] = temp;
      temp = temp/(9*mads[j]);
      rUs[i + rn*j] = temp;
      rAs[i + rn*j] = 1;
      if(temp < -1 | temp > 1)
        rAs[i + rn*j] = 0;
    };
  };
  
  for(j=0; j < cn; j++)
  {
    tot1 = 0;
    for(i=0; i < rn; i++)
    {
      s1 = rAs[i + rn*j]*rXMs[i + rn*j]*(1-rUs[i + rn*j]*rUs[i + rn*j])*(1-rUs[i + rn*j]*rUs[i + rn*j]);
      tot1 += s1*s1;
    }
    top = rn*tot1;
    
    tot1 = 0;
    tot2 = 0;
    for(i=0; i < rn; i++)
    {
      s1 = rAs[i + rn*j]*(1-rUs[i + rn*j]*rUs[i + rn*j])*(1-5*rUs[i + rn*j]*rUs[i + rn*j]);
      tot1 += s1;
      tot2 += s1;
    }
    bot = tot1*tot2;
    
    rbidiag[j] = top/bot;
  }
  
  int ii, jj;
  
  for(ii=0; ii < (cn-1); ii++)
  {
    for(jj=(ii+1); jj < cn; jj++)
    {
      tot1 = 0;
      for(i=0; i < rn; i++)
      {
        s1 = rAs[i + rn*ii]*rXMs[i + rn*ii]*(1-rUs[i + rn*ii]*rUs[i + rn*ii])*(1-rUs[i + rn*ii]*rUs[i + rn*ii]);
        s2 = rAs[i + rn*jj]*rXMs[i + rn*jj]*(1-rUs[i + rn*jj]*rUs[i + rn*jj])*(1-rUs[i + rn*jj]*rUs[i + rn*jj]);
        tot1 += s1*s2;
      }
      top = rn*tot1;
  
      tot1 = 0;
      tot2 = 0;
      for(i=0; i < rn; i++)
      {
        s1 = rAs[i + rn*ii]*(1-rUs[i + rn*ii]*rUs[i + rn*ii])*(1-5*rUs[i + rn*ii]*rUs[i + rn*ii]);
        s2 = rAs[i + rn*jj]*(1-rUs[i + rn*jj]*rUs[i + rn*jj])*(1-5*rUs[i + rn*jj]*rUs[i + rn*jj]);
        tot1 += s1;
        tot2 += s2;
      }
      bot = tot1*tot2;
  
      rbicor[ii + cn*jj] = (top/bot)/sqrt(rbidiag[ii]*rbidiag[jj]);
      rbicor[jj + cn*ii] = rbicor[ii + cn*jj];
    }
  }
  
  for(j=0; j < cn; j++)
    rbicor[j + cn*j] = 1;
  
  UNPROTECT(5);
  return(bicor);
}

R_CallMethodDef callMethods[]  = {
    {"f0worker", (DL_FUNC)&f0worker, 4},
    {"bwmcCworker", (DL_FUNC)&bwmcCworker, 5},
    {NULL, NULL, 0}
};

void R_init_EBcoexpress(DllInfo *info)
{
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
