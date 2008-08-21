#include "ljrlib.h"
void test12(double *y,double *n,double *tm,double *x,double *ofst,int *Nptr,int *mptr,int *Rptr,double *p){
 int N=*Nptr;
 int m=*mptr;
 int R=*Rptr;

 double *b1;
 double *b;
 double *simy;

 double g1[2], g[3];
 double tau1[1], tau[2];
 double liko1,liko2,lik1,lik2;
 double obs,sim;
 int count=0;
 int i;

 b1=malloc((m+1)*sizeof(double));
 b=malloc((m+1)*sizeof(double));
 simy=malloc(N*sizeof(double));

 ljr1(y,n,tm,x,ofst,b1,g1,tau1,&N,&m,&liko1);
 ljr2(y,n,tm,x,ofst,b,g,tau,&N,&m,&liko2);

 count=0;
 obs=liko2-liko1;
 for (i=0;i<=R;i++){
  rgy(b1,g1,tau1,n,tm,x,ofst,simy,N,m,1);
  ljr1(simy,n,tm,x,ofst,b,g,tau,&N,&m,&lik1);
  ljr2(simy,n,tm,x,ofst,b,g,tau,&N,&m,&lik2);
  sim=lik2-lik1;
  if (sim>obs)
   count++;
 }
 p[0]=count/(R+0.0);
 free(b1);
 free(b);
 free(simy); 
}
