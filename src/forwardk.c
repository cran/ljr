#include "ljrlib.h"
void forwardk(double *y,double *n,double *tm,double *x,double *ofst,int *Nptr,int *mptr,int *Rptr,double *p,double *alphaptr,int *ncps){
 int N=*Nptr;
 int m=*mptr;
 int R=*Rptr;
 double alpha=*alphaptr;

 double *b0;
 double *b;
 double *simy;
 double *g0;
 double *g;
 double *tau0;
 double *tau;

 double obs,sim;
 int count;
 int i,j;
 double lik0,lik1,liko0,liko1;
 int nullncps,altncps;

 b0=malloc((m+1)*sizeof(double));
 b=malloc((m+1)*sizeof(double));
 simy=malloc(N*sizeof(double));
 g0=malloc((N-2)*sizeof(double));
 g=malloc((N-1)*sizeof(double));
 tau0=malloc((N-3)*sizeof(double));
 tau=malloc((N-2)*sizeof(double));

 nullncps=0;
 altncps=1;
 p[0]=-1;
 j=0;
 while (p[j]<alpha){
  if (nullncps==0)
   ljr0(y,n,tm,x,ofst,b0,g0,&N,&m,&liko0);
  else
   ljrk(&nullncps,y,n,tm,x,ofst,b0,g0,tau0,&N,&m,&liko0);
  ljrk(&altncps,y,n,tm,x,ofst,b,g,tau,&N,&m,&liko1);
  obs=liko1-liko0;
  count=0;
  for (i=0;i<R;i++){
   rgy(b0,g0,tau,n,tm,x,ofst,simy,N,m,nullncps);
   if (nullncps==0)
    ljr0(simy,n,tm,x,ofst,b,g,&N,&m,&lik0);
   else
    ljrk(&nullncps,simy,n,tm,x,ofst,b,g,tau,&N,&m,&lik0);
   ljrk(&altncps,simy,n,tm,x,ofst,b,g,tau,&N,&m,&lik1);
   sim=lik1-lik0;
   if (sim>obs)
    count++;
  }
  p[j]=count/(R+0.0);
  if (p[j]<alpha){
   altncps++;
   nullncps++;
   j++;
   p[j]=-1;
  }
 }
 *ncps=nullncps;

 free(b0);
 free(b);
 free(simy);
 free(g0);
 free(g);
 free(tau0);
 free(tau);
}
