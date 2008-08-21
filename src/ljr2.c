#include "ljrlib.h"
void ljr2(double *y,double *n,double *tm,double *x,double *ofst,double *beta,double *gamma,double *tau,int *Nptr,int *mptr,double *zlik){
 int N=*Nptr;
 int m=*mptr;

 int i,i1,i2;
 int ncps[2]={2,0};

 double *z;
 int *dt;
 double *b; 

 int ifixed[2]={N,N}; 
 int ifree[2]={N,N};
 int ndt[1];
 double zliktemp[1];
 double temptau,temptau2;
 int tempintptr;

 z=malloc(N*(m+6)*sizeof(double));
 dt=malloc(N*sizeof(int));
 b=malloc((m+6)*sizeof(double));

 for (i=0;i<m+6;i++)
  b[i]=0;

 rgi(tm,dt,ndt,N);	
 ifixed[0]=dt[2];
 ifixed[1]=dt[3];
 rz(z,x,tm,ifixed,ifree,ncps,N,m,1);
 tempintptr=m+4;
 lr(y,n,z,ofst,b,N,tempintptr,zlik);
 rcf(beta,gamma,b,m,ncps);
 tau[0]=tm[dt[2]-1];
 tau[1]=tm[dt[3]-1];
 for (i2=4;i2<*ndt;i2++){
  mvfixed(z,tm,ncps,dt,N,m,2,i2);
  lr(y,n,z,ofst,b,N,tempintptr,zliktemp);
  if (*zliktemp>*zlik){
   rcf(beta,gamma,b,m,ncps);
   tau[1]=tm[dt[i2]-1];
   *zlik=*zliktemp;
  }
 }
 for (i1=3;i1<*ndt-1;i1++){
  mvfixed(z,tm,ncps,dt,N,m,1,i1);
  fillfixed(z,tm,ncps,dt,N,m,2,i1+1,i1+1,*ndt);
  lr(y,n,z,ofst,b,N,tempintptr,zliktemp);
  if (*zliktemp>*zlik){
   rcf(beta,gamma,b,m,ncps);
   tau[0]=tm[dt[i1]-1];
   tau[1]=tm[dt[i2]-1];
   *zlik=*zliktemp;
  }
  for (i2=i1+2;i2<*ndt;i2++){
   mvfixed(z,tm,ncps,dt,N,m,2,i2);
   lr(y,n,z,ofst,b,N,tempintptr,zliktemp);
   if (*zliktemp>*zlik){
    rcf(beta,gamma,b,m,ncps);
    tau[0]=tm[dt[i1]-1];
    tau[1]=tm[dt[i2]-1];
    *zlik=*zliktemp;
   }
  }
 }

 ncps[0]=1;
 ncps[1]=1;
 tempintptr=m+5;
 ifixed[0]=dt[2];
 ifree[0]=dt[3];
 rz(z,x,tm,ifixed,ifree,ncps,N,m,0);
 lr(y,n,z,ofst,b,N,tempintptr,zliktemp);
 if (*zliktemp>*zlik){
  if (b[m+1]!=b[m+3])
   temptau=(b[m+2]-b[m])/(b[m+1]-b[m+3]);
  else
   temptau=tm[0];
  if ((temptau>tm[dt[3]-1])&(temptau<tm[dt[3]])){
   rcf(beta,gamma,b,m,ncps);
   tau[0]=temptau;
   tau[1]=tm[dt[2]-1];
   *zlik=*zliktemp;
  }
 }    
 for (i2=4;i2<*ndt-1;i2++){
  mvfree(z,tm,dt,N,m,1,i2);
 lr(y,n,z,ofst,b,N,tempintptr,zliktemp);
 if (*zliktemp>*zlik){
  if (b[m+1]!=b[m+3])
   temptau=(b[m+2]-b[m])/(b[m+1]-b[m+3]);
  else
   temptau=tm[0];
  if ((temptau>tm[dt[i2]-1])&(temptau<tm[dt[i2]])){
   rcf(beta,gamma,b,m,ncps);
   tau[0]=temptau;
   tau[1]=tm[dt[2]-1];
   *zlik=*zliktemp;
  }
 }    
 }
 for (i1=3;i1<*ndt;i1++){
  mvfixed(z,tm,ncps,dt,N,m,1,i1);
  ifree[0]=dt[2];
  fillfree(z,tm,ncps,dt,N,m,ifree);
  if (i1>3){
   lr(y,n,z,ofst,b,N,tempintptr,zliktemp);
   if (*zliktemp>*zlik){
    if (b[m+1]!=b[m+3])
     temptau=(b[m+2]-b[m])/(b[m+1]-b[m+3]);
    else
     temptau=tm[0];
    if ((temptau>tm[dt[2]-1])&(temptau<tm[dt[2]])){
     rcf(beta,gamma,b,m,ncps);
     tau[0]=temptau;
     tau[1]=tm[dt[i1]-1];
     *zlik=*zliktemp;
    }
   }  
  }
  for (i2=3;i2<*ndt-1;i2++){
   mvfree(z,tm,dt,N,m,1,i2);
   if ((i2<i1-1)|(i2>i1)){
    lr(y,n,z,ofst,b,N,tempintptr,zliktemp);
    if (*zliktemp>*zlik){
     if (b[m+1]!=b[m+3])
      temptau=(b[m+2]-b[m])/(b[m+1]-b[m+3]);
     else
      temptau=tm[0];
     if ((temptau>tm[dt[i2]-1])&(temptau<tm[dt[i2]])){
      rcf(beta,gamma,b,m,ncps);
      tau[0]=temptau;
      tau[1]=tm[dt[i1]-1];
      *zlik=*zliktemp;
     }
    }
   }
  }
 }
 
 ncps[0]=0;
 ncps[1]=2;
 tempintptr=m+6;
 for (i1=2;i1<*ndt-3;i1++){
  ifree[0]=dt[i1];
  ifree[1]=dt[i1+2];
  rz(z,x,tm,ifixed,ifree,ncps,N,m,0);
  lr(y,n,z,ofst,b,N,tempintptr,zliktemp);
  if (*zliktemp>*zlik){
   if (b[m+1]!=b[m+3])
    temptau=(b[m+2]-b[m])/(b[m+1]-b[m+3]);
   else
    temptau=tm[0];
   if (b[m+3]!=b[m+5])
    temptau2=(b[m+4]-b[m+2])/(b[m+3]-b[m+5]);
   else
    temptau2=tm[0];
   if ((temptau>tm[dt[i1]-1])&(temptau<tm[dt[i1]])&(temptau2>tm[dt[i1+2]-1])&(temptau2<tm[dt[i1+2]])){
    rcf(beta,gamma,b,m,ncps);
    tau[0]=temptau;
    tau[1]=temptau2;
    *zlik=*zliktemp;
   }
  }
  for (i2=i1+3;i2<*ndt-1;i2++){
   mvfree(z,tm,dt,N,m,2,i2);
   lr(y,n,z,ofst,b,N,tempintptr,zliktemp);
   if (*zliktemp>*zlik){
    if (b[m+1]!=b[m+3])
     temptau=(b[m+2]-b[m])/(b[m+1]-b[m+3]);
    else
     temptau=tm[0];
    if (b[m+3]!=b[m+5])
     temptau2=(b[m+4]-b[m+2])/(b[m+3]-b[m+5]);
    else
     temptau2=tm[0];
    if ((temptau>tm[dt[i1]-1])&(temptau<tm[dt[i1]])&(temptau2>tm[dt[i2]-1])&(temptau2<tm[dt[i2]])){
     rcf(beta,gamma,b,m,ncps);
     tau[0]=temptau;
     tau[1]=temptau2;
     *zlik=*zliktemp;
    }
   }
  }
 }
 free(z);
 free(dt);
 free(b);  
}
