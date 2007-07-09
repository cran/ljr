#include "ljrlib.h"
void rgy(double *beta,double *gamma,double *tau,double *n,double *tm,double *x,double *ofst,double *y,int N,int m,int ncps){
 int i,j;
 int nint;
 double eta,p;
 double rb[1];
 for (i=0;i<N;i++){
  eta=ofst[i]+beta[0]+gamma[0]*tm[i];
  for (j=0;j<m;j++)
   eta+=beta[j+1]*x[i+j*N];
  for (j=0;j<ncps;j++)
   eta+=gamma[j+1]*fmax2(tm[i]-tau[j],0);
  p=1-1/(1+exp(eta));
  nint=(int)floor(n[i]+.1);
  rrbinom(&nint,&p,rb);
  y[i]=rb[0];
 }
}
