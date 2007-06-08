#include "ljrlib.h"
void test01(double *y,double *n,double *tm,double *x,double *ofst,int *Nptr,int *mptr,int *Rptr,double *p){
 int N=*Nptr;
 int m=*mptr;
 int R=*Rptr;

 double *b0;
 double *b;
 double *simy;

 double g0[1], g[2];
 double tau[1];
 double lik0,lik1;
 double obslambda,simlambda;
 int count=0;
 int i;

 b0=malloc((m+1)*sizeof(double));
 b=malloc((m+1)*sizeof(double));
 simy=malloc(N*sizeof(double));
// estimate parameters for null model and computes observed LRT
 ljr0(y,n,tm,x,ofst,b0,g0,&N,&m,&lik0);
 ljr1(y,n,tm,x,ofst,b,g,tau,&N,&m,&lik1);
 obslambda=lik1-lik0;
// estimates p-value via Monte Carlo
 for (i=0;i<R;i++){
  rgy(b0,g0,tau,n,tm,x,ofst,simy,N,m,0);
  ljr0(simy,n,tm,x,ofst,b,g,&N,&m,&lik0);
  ljr1(simy,n,tm,x,ofst,b,g,tau,&N,&m,&lik1);
  simlambda=lik1-lik0;
  if (simlambda>obslambda)
   count++;
 }
 p[0]=count/(R+0.0);
 free(b0);
 free(b);
 free(simy);
}
