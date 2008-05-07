#include "ljrlib.h"
void ljr1rmint(double *y,double *n,double *tm,double *x,double *ofst,double *beta,double *gamma,double *tau,int *Nptr,int *mptr,double *zlik){
 int N=*Nptr;
 int m=*mptr;

 int i;
 int ncps[2]={1,0};

 double *z;
 int *dt;
 double *b;

 int ifixed[1]={N};
 int ifree[1]={N};
 int ndt[1];
 double zliktemp[1];
 double temptau;
 int tempintptr;

 z=malloc(N*(m+4)*sizeof(double));
 dt=malloc(N*sizeof(int));
 b=malloc((m+4)*sizeof(double));
 
 for (i=0;i<m+4;i++)
  b[i]=0;

 rgi(tm,dt,ndt,N);
 ifixed[0]=dt[2];
 rz(z,x,tm,ifixed,ifree,ncps,N,m,1);
 rzrmrow(z,m+1,N,m+3);
 tempintptr=m+2;
 lr(y,n,z,ofst,b,N,tempintptr,zlik);
 rcfrmint(beta,gamma,b,m,ncps);
 tau[0]=tm[dt[2]-1];

 for (i=3;i<*ndt;i++){
  mvfixedrm(z,tm,ncps,dt,N,m,1,i);
  lr(y,n,z,ofst,b,N,tempintptr,zliktemp);
  if (*zliktemp>*zlik){
   rcfrmint(beta,gamma,b,m,ncps);
   tau[0]=tm[dt[i]-1];
   *zlik=*zliktemp;
  }
 }
 
 ncps[0]=0;
 ncps[1]=1;
 ifree[0]=dt[2];
 rz(z,x,tm,ifixed,ifree,ncps,N,m,0);
 rzrmrow(z,m+1,N,m+4);
 tempintptr=m+3;
 lr(y,n,z,ofst,b,N,tempintptr,zliktemp);
 if (*zliktemp>*zlik){
  if (b[m]!=b[m+2])
   temptau=b[m+1]/(b[m]-b[m+2]);
  else
   temptau=tm[0];
  if ((temptau>tm[dt[2]-1])&(temptau<tm[dt[2]])){
   rcfrmint(beta,gamma,b,m,ncps);
   tau[0]=temptau;
   *zlik=*zliktemp;
  }  
 }

 for (i=3;i<*ndt-1;i++){
  mvfreermint(z,tm,dt,N,m,1,i);
  lr(y,n,z,ofst,b,N,tempintptr,zliktemp);
  if (*zliktemp>*zlik){
   if (b[m]!=b[m+2])
    temptau=b[m+1]/(b[m]-b[m+2]);
   else
    temptau=tm[0];
   if ((temptau>tm[dt[i]-1])&(temptau<tm[dt[i]])){
    rcfrmint(beta,gamma,b,m,ncps);
    tau[0]=temptau;
    *zlik=*zliktemp;
   }  
  }
 }
 free(z);
 free(dt);
 free(b);
}
