#include "ljrlib.h"
void ljrkk(int *kptr,double *y,double *n,double *tm,double *x,double *ofst,int *Nptr,int *mptr,int *Rptr,double *p){
 int k=*kptr;
 int N=*Nptr;
 int m=*mptr;
 int R=*Rptr;

 double *b0;
 double *b;
 double *g0;
 double *g;
 double *tau0;
 double *tau; 
 double *simy;
 double *xx;

 double liko0,liko1,lik0,lik1;
 double obs,sim;
 int count;
 int i,j,h;
 int mm=m-1;

 b0=malloc((m+1)*sizeof(double));
 b=malloc((m+1)*sizeof(double));
 g0=malloc((k+1)*sizeof(double));
 g=malloc((k+1)*sizeof(double));
 tau0=malloc(k*sizeof(double));
 tau=malloc(k*sizeof(double));
 simy=malloc(N*sizeof(double));
 xx=malloc(N*m*sizeof(double));

 for (i=0;i<m;i++){
  count=0;
  for (h=0;h<m*N;h++)
   xx[h]=x[h];
  rzrmrow(xx,i+1,N,m);
  ljrk(&k,y,n,tm,xx,ofst,b0,g0,tau0,&N,&mm,&liko0);
  ljrk(&k,y,n,tm,x,ofst,b,g,tau,&N,&m,&liko1);
  obs=liko1-liko0;
  for (j=0;j<R;j++){
   rgy(b0,g0,tau0,n,tm,xx,ofst,simy,N,mm,k);
   ljrk(&k,simy,n,tm,xx,ofst,b,g,tau,&N,&mm,&lik0);
   ljrk(&k,simy,n,tm,x,ofst,b,g,tau,&N,&m,&lik1);
   sim=lik1-lik0;
   if (sim>obs)
    count++;
  }
  p[i]=count/(R+0.0);
 }
 count=0;
//iprpm('k',&k,1,1);
//prpm('y',y,N,1);
//prpm('n',n,N,1);
//prpm('t',tm,N,1);
//prpm('o',ofst,N,1);
 ljrkrmint(&k,y,n,tm,x,ofst,b0,g0,tau0,&N,&m,&liko0);
//prpm('b',b0,m+1,1);
//prpm('g',g0,k+1,1);
//prpm('l',&liko0,1,1);

 ljrk(&k,y,n,tm,x,ofst,b,g,tau,&N,&m,&liko1);
 obs=liko1-liko0;
 for (j=0;j<R;j++){
  rgy(b0,g0,tau0,n,tm,x,ofst,simy,N,m,k);
  ljrkrmint(&k,simy,n,tm,x,ofst,b,g,tau,&N,&m,&lik0);
  ljrk(&k,simy,n,tm,x,ofst,b,g,tau,&N,&m,&lik1);
  sim=lik1-lik0;
  if (sim>obs)
   count++;
 }
 p[m]=count/(R+0.0);
 count=0;
 ljrkrmtm(&k,y,n,tm,x,ofst,b0,g0,tau0,&N,&m,&liko0);
 ljrk(&k,y,n,tm,x,ofst,b,g,tau,&N,&m,&liko1);
 obs=liko1-liko0;
 for (j=0;j<R;j++){
  rgy(b0,g0,tau0,n,tm,x,ofst,simy,N,m,k);
  ljrkrmtm(&k,simy,n,tm,x,ofst,b,g,tau,&N,&m,&lik0);
  ljrk(&k,simy,n,tm,x,ofst,b,g,tau,&N,&m,&lik1);
  sim=lik1-lik0;
  if (sim>obs)
   count++;
 }
 p[m+1]=count/(R+0.0);
 free(b0);
 free(b);
 free(g0);
 free(g);
 free(tau0);
 free(tau);
 free(simy);
 free(xx);
}
