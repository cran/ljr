#include "ljrlib.h"
void testjk(int *jptr,int *kptr,double *y,double *n,double *tm,double *x,double *ofst,int *Nptr,int *mptr,int *Rptr,double *p){
 int j=*jptr;
 int k=*kptr;
 int N=*Nptr;
 int m=*mptr;
 int R=*Rptr;
 int i;

 if (j>k){
  i=j;
  j=k;
  k=i;
 }

 double *b0;
 double *b;
 double *simy;
 double *g0;
 double *g;
 double *tau0;
 double *tau;

 double obs,sim;
 int count;
 double lik0,lik1,liko0,liko1;

 b0=malloc((m+1)*sizeof(double));
 b=malloc((m+1)*sizeof(double));
 simy=malloc(N*sizeof(double));
 g0=malloc((j+1)*sizeof(double));
 g=malloc((k+1)*sizeof(double));
 tau0=malloc(j*sizeof(double));
 tau=malloc(k*sizeof(double));

 ljrk(&j,y,n,tm,x,ofst,b0,g0,tau0,&N,&m,&liko0);
 ljrk(&k,y,n,tm,x,ofst,b,g,tau,&N,&m,&liko1);
 obs=liko1-liko0;

 count=0;
 for (i=0;i<=R;i++){
  rgy(b0,g0,tau0,n,tm,x,ofst,simy,N,m,j);
  ljrk(&j,simy,n,tm,x,ofst,b,g,tau,&N,&m,&lik0);
  ljrk(&k,simy,n,tm,x,ofst,b,g,tau,&N,&m,&lik1);
  sim=lik1-lik0;
  if (sim>obs)
   count++;
 }
 p[0]=count/(R+0.0);

 free(b0);
 free(b);
 free(simy);
 free(g0);
 free(g);
 free(tau0);
 free(tau);
}
