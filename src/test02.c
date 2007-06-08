#include "ljrlib.h"
void test02(double *y,double *n,double *tm,double *x,double *ofst,int *Nptr,int *mptr,int *Rptr,double *p){
 int N=*Nptr;
 int m=*mptr;
 int R=*Rptr;

 double *b0;
 double *b;
 double *simy;

 double g0[1], g[3];
 double tau[2];
 double liko0,liko2,lik0,lik2;
 double obs,sim;
 int count=0;
 int i;

 b0=malloc((m+1)*sizeof(double));
 b=malloc((m+1)*sizeof(double));
 simy=malloc(N*sizeof(double));

 ljr0(y,n,tm,x,ofst,b0,g0,&N,&m,&liko0);
 ljr2(y,n,tm,x,ofst,b,g,tau,&N,&m,&liko2);
 obs=liko2-liko0;
 for (i=0;i<=R;i++){
  rgy(b0,g0,tau,n,tm,x,ofst,simy,N,m,0);
  ljr0(simy,n,tm,x,ofst,b,g,&N,&m,&lik0);
  ljr2(simy,n,tm,x,ofst,b,g,tau,&N,&m,&lik2);
  sim=lik2-lik0;
  if (sim>obs)
   count++;
 }
 p[0]=count/(R+0.0);
 free(b0);
 free(b);
 free(simy); 
}
