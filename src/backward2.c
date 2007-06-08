#include "ljrlib.h"
void backward2(double *y,double *n,double *tm,double *x,double *ofst,int *Nptr,int *mptr,int *Rptr,double *p,double *alphaptr,int *nulls,int *alts){
 int N=*Nptr;
 int m=*mptr;
 int R=*Rptr;
 double alpha=*alphaptr;

 double *b0;
 double *b1;
 double *b;
 double *simy;

 double g0[1];
 double g1[2];
 double g[3];
 double tau1[1];
 double tau[2];
	 
 double liko0,liko1,liko2,lik0,lik1,lik2;
 double obs,sim;
 int count=0;
 int i;
 
 b0=malloc((m+1)*sizeof(double));
 b1=malloc((m+1)*sizeof(double));
 b=malloc((m+1)*sizeof(double));
 simy=malloc(N*sizeof(double));

// estimate parameters for null model and computes observed LRTs
 ljr0(y,n,tm,x,ofst,b0,g0,&N,&m,&liko0);
 ljr1(y,n,tm,x,ofst,b1,g1,tau1,&N,&m,&liko1);
 ljr2(y,n,tm,x,ofst,b,g,tau,&N,&m,&liko2);
// Test 0 vs 2.
 nulls[0]=0;
 alts[0]=2;
 obs=liko2-liko0;
 for (i=0;i<R;i++){
  rgy(b0,g0,tau,n,tm,x,ofst,simy,N,m,0);
  ljr0(simy,n,tm,x,ofst,b,g,&N,&m,&lik0);
  ljr2(simy,n,tm,x,ofst,b,g,tau,&N,&m,&lik2);
  sim=lik2-lik0;
  if (sim>obs)
   count++;
 }
 p[0]=count/(R+0.0);
 if (p[0]>alpha/2){
// If 0 vs 2 is not rejected, test 0 vs 1.
  count=0;
  nulls[1]=0;
  alts[1]=1;
  obs=liko1-liko0;
  for (i=0;i<R;i++){
   rgy(b0,g0,tau,n,tm,x,ofst,simy,N,m,0);
   ljr0(simy,n,tm,x,ofst,b,g,&N,&m,&lik0);
   ljr1(simy,n,tm,x,ofst,b,g,tau,&N,&m,&lik1);
   sim=lik1-lik0;
   if (sim>obs)
    count++;
  }
  p[1]=count/(R+0.0);
// If 0 vs 1 is not rejected, then use 0 joinpoints.  Otherwise use 1.
  if (p[1]>alpha/2)
   nulls[2]=0;
  else
   nulls[2]=1;
 }
 else{
// If 0 vs 2 is rejected, test 1 vs 2.
  count=0;
  nulls[1]=1;
  alts[1]=2;
  obs=liko2-liko1;
  for (i=0;i<R;i++){
   rgy(b1,g1,tau1,n,tm,x,ofst,simy,N,m,1);
   ljr1(simy,n,tm,x,ofst,b,g,tau,&N,&m,&lik1);
   ljr2(simy,n,tm,x,ofst,b,g,tau,&N,&m,&lik2);
   sim=lik2-lik1;
   if (sim>obs)
    count++;
  }
  p[1]=count/(R+0.0);
// If 1 vs 2 is not rejected, then use 1 joinpoint.  Otherwise use 2.
  if (p[1]>alpha/2)
   nulls[2]=1;
  else
   nulls[2]=2;
 }
 free(b0);
 free(b1);
 free(b);
 free(simy);
}
