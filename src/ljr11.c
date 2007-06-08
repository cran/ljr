#include "ljrlib.h"
void ljr11(double *y,double *n,double *tm,double *x,double *ofst,int *Nptr,int *mptr,int *Rptr,double *p){
 int N=*Nptr;
 int m=*mptr;
 int R=*Rptr;

 double *b1;
 double *b;
 double *simy;
 double *xx;

 double g1[2], g[2];
 double tau1[1], tau[1];
 double liko0,liko1,lik0,lik1;
 double obs,sim;
 int count;
 int i,j,k;
 int mm=m-1;

 b1=malloc((m+1)*sizeof(double));
 b=malloc((m+1)*sizeof(double));
 simy=malloc(N*sizeof(double));
 xx=malloc(N*m*sizeof(double));

// estimate parameters for null model and compute observed LRTs
 for (i=0;i<m;i++){
  count=0;
  for (k=0;k<m*N;k++)
   xx[k]=x[k];
  rzrmrow(xx,i+1,N,m);
  ljr1(y,n,tm,xx,ofst,b1,g1,tau1,&N,&mm,&liko0);
  ljr1(y,n,tm,x,ofst,b,g,tau,&N,&m,&liko1);
  obs=liko1-liko0;
  for (j=0;j<R;j++){
   rgy(b1,g1,tau1,n,tm,xx,ofst,simy,N,mm,1);
   ljr1(simy,n,tm,xx,ofst,b,g,tau,&N,&mm,&lik0);
   ljr1(simy,n,tm,x,ofst,b,g,tau,&N,&m,&lik1);
   sim=lik1-lik0;
   if (sim>obs)
    count++;
  }
  p[i]=count/(R+0.0);
 }
 count=0;
 ljr1rmint(y,n,tm,x,ofst,b1,g1,tau1,&N,&m,&liko0);
 ljr1(y,n,tm,x,ofst,b,g,tau,&N,&m,&liko1);
 obs=liko1-liko0;
 for (j=0;j<R;j++){
  rgy(b1,g1,tau1,n,tm,x,ofst,simy,N,m,1);
  ljr1rmint(simy,n,tm,x,ofst,b,g,tau,&N,&m,&lik0);
  ljr1(simy,n,tm,x,ofst,b,g,tau,&N,&m,&lik1);
  sim=lik1-lik0;
  if (sim>obs)
   count++;
 }
 p[m]=count/(R+0.0);
 count=0;
 ljr1rmtm(y,n,tm,x,ofst,b1,g1,tau1,&N,&m,&liko0);
 ljr1(y,n,tm,x,ofst,b,g,tau,&N,&m,&liko1);
 obs=liko1-liko0;
 for (j=0;j<R;j++){
  rgy(b1,g1,tau1,n,tm,x,ofst,simy,N,m,1);
  ljr1rmtm(simy,n,tm,x,ofst,b,g,tau,&N,&m,&lik0);
  ljr1(simy,n,tm,x,ofst,b,g,tau,&N,&m,&lik1);
  sim=lik1-lik0;
  if (sim>obs)
   count++;
 }
 p[m+1]=count/(R+0.0);
 free(b1);
 free(b);
 free(simy);
 free(xx);
}
