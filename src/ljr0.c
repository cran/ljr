#include "ljrlib.h"
void ljr0(double *y,double *n,double *tm,double *x,double *ofst,double *beta,double *gamma,int *Nptr,int *mptr,double *zlik){
 int N=*Nptr;
 int m=*mptr;
 int i;

 double *z;
 double *b;
 int ifixed[1]={0};
 int ifree[1]={0};
 int ncps[2]={0,0};
 int tempintptr;

 z=malloc(N*(m+2)*sizeof(double));
 b=malloc((m+2)*sizeof(double));
 rz(z,x,tm,ifixed,ifree,ncps,N,m,1);
 for (i=0;i<m+2;i++)
  b[i]=0;
 tempintptr=m+2;
 lr(y,n,z,ofst,b,N,tempintptr,zlik);
 rcf(beta,gamma,b,m,ncps);
 free(z);
 free(b);
}
