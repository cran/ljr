#include "ljrlib.h"
void ljr22(double *y,double *n,double *tm,double *x,double *ofst,int *Nptr,int *mptr,int *Rptr,double *p){
 int N=*Nptr;
 int m=*mptr;
 int R=*Rptr;

 double *b2;
 double *b;
 double *simy;
 double *xx;

 double g2[3], g[3];
 double tau2[2], tau[2];
 double liko0,liko1,lik0,lik1;
 double obs,sim;
 int count;
 int i,j,k;
 int mm=m-1;

 b2=malloc((m+1)*sizeof(double));
 b=malloc((m+1)*sizeof(double));
 simy=malloc(N*sizeof(double));
 xx=malloc(N*m*sizeof(double));

 for (i=0;i<m;i++){
  count=0;
  for (k=0;k<m*N;k++)
   xx[k]=x[k];
  rzrmrow(xx,i+1,N,m);
  ljr2(y,n,tm,xx,ofst,b2,g2,tau2,&N,&mm,&liko0);
  ljr2(y,n,tm,x,ofst,b,g,tau,&N,&m,&liko1);
  obs=liko1-liko0;
  for (j=0;j<R;j++){
   rgy(b2,g2,tau2,n,tm,xx,ofst,simy,N,mm,2);
   ljr2(simy,n,tm,xx,ofst,b,g,tau,&N,&mm,&lik0);
   ljr2(simy,n,tm,x,ofst,b,g,tau,&N,&m,&lik1);
   sim=lik1-lik0;
   if (sim>obs)
    count++;
  }
  p[i]=count/(R+0.0);
 }
 count=0;
 ljr2rmint(y,n,tm,x,ofst,b2,g2,tau2,&N,&m,&liko0);
 ljr2(y,n,tm,x,ofst,b,g,tau,&N,&m,&liko1);
 obs=liko1-liko0;
 for (j=0;j<R;j++){
  rgy(b2,g2,tau2,n,tm,x,ofst,simy,N,m,2);
  ljr2rmint(simy,n,tm,x,ofst,b,g,tau,&N,&m,&lik0);
  ljr2(simy,n,tm,x,ofst,b,g,tau,&N,&m,&lik1);
  sim=lik1-lik0;
  if (sim>obs)
   count++;
 }
 p[m]=count/(R+0.0);
 count=0;
 ljr2rmtm(y,n,tm,x,ofst,b2,g2,tau2,&N,&m,&liko0);
 ljr2(y,n,tm,x,ofst,b,g,tau,&N,&m,&liko1);
 obs=liko1-liko0;
 for (j=0;j<R;j++){
  rgy(b2,g2,tau2,n,tm,x,ofst,simy,N,m,2);
  ljr2rmtm(simy,n,tm,x,ofst,b,g,tau,&N,&m,&lik0);
  ljr2(simy,n,tm,x,ofst,b,g,tau,&N,&m,&lik1);
  sim=lik1-lik0;
  if (sim>obs)
   count++;
 }
 p[m+1]=count/(R+0.0);
 free(b2);
 free(b);
 free(simy);
 free(xx);
}

