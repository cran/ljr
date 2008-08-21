#include "ljrlib.h"
void ljrkrmtm(int *kptr,double *y,double *n,double *tm,double *x,double *ofst,double *beta,double *gamma,double *tau,int *Nptr,int *mptr,double *zlik){
 int k=*kptr;
 int N=*Nptr;
 int m=*mptr;

 int i;
 int ncps[2]={k,0};
 int ndt[1];
 double zliktemp[1];
 int tempintptr;
 int M;
 int v;
 int rst=0;
 int good;
 int numfixed[1];
 int nojpinint;

 double *z;
 int *ifixed;
 int *ifree;
 int *dt;
 double *b;
 double *temptau;
 int *ik;

 z=malloc(N*(m+2*k+2)*sizeof(double));
 ifixed=malloc(k*sizeof(int));
 ifree=malloc(k*sizeof(int));
 dt=malloc(N*sizeof(int));
 b=malloc((m+2*k+2)*sizeof(double));
 temptau=malloc(k*sizeof(double));
 ik=malloc(k*sizeof(int));

  for (i=0;i<m+2*k+2;i++)
   b[i]=0;

  rgi(tm,dt,ndt,N);
  M=*ndt-2;
  if (k>M)
   k=M;
  v=k-1;
  *numfixed=k;
  for (i=0;i<k;i++)
   ik[i]=1+i;
  srz(ik,dt,ifixed,ifree,ncps,k,M,*numfixed);
  rz(z,x,tm,ifixed,ifree,ncps,N,m,1);
  rzrmrow(z,m+2,N,m+k+2+ncps[1]);
  tempintptr=m+k+1+ncps[1];
  lr(y,n,z,ofst,b,N,tempintptr,zlik);
  rcfrmtm(beta,gamma,b,m,ncps);
  for (i=0;i<k;i++)
   tau[i]=tm[ifixed[i]-1];
  while (v>-1){
   if (ik[v]<2*M+v-k){
    if (rst==1){
     ik[v]++;
     for (i=v+1;i<k;i++)
      ik[i]=ik[i-1]+1;
     v=k-1;
     rst=0;
    }
    else
     ik[v]++;
    good=checkseq(ik,M,k,numfixed);
    if (good){
     srz(ik,dt,ifixed,ifree,ncps,k,M,*numfixed);
     rz(z,x,tm,ifixed,ifree,ncps,N,m,0);
     rzrmrow(z,m+2,N,m+k+2+ncps[1]);
     tempintptr=m+k+1+ncps[1];
     lr(y,n,z,ofst,b,N,tempintptr,zliktemp);
     if (*zliktemp>*zlik){
      if (b[m+2]!=0){
       temptau[0]=(b[m]-b[m+1])/b[m+2];
       nojpinint=0;
      }
      else
       nojpinint=1;
      nojpinint=0;
      i=1;
      while ((nojpinint==0)&(i<ncps[1])){
       if (b[m+2*i]!=b[m+2+2*i])
        temptau[i]=(b[m+1+2*i]-b[m-1+2*i])/(b[m+2*i]-b[m+2+2*i]);
       else
        nojpinint=1;
       i++;
      }
      i=0;
      while ((nojpinint==0)&(i<ncps[1])){
       if ((temptau[i]<=tm[ifree[i]-1])|(temptau[i]>=tm[ifree[i]]))
        nojpinint=1;
       i++;
      }
      if (nojpinint==0){
       rcfrmtm(beta,gamma,b,m,ncps);
       for (i=0;i<ncps[0];i++)
        tau[i]=tm[ifixed[i]-1];
       for (i=0;i<k-ncps[0];i++)
        tau[ncps[0]+i]=temptau[i];
       *zlik=*zliktemp;
      }
     }
    }
   }
   else{
    v--;
    rst=1;
   }
  }
 free(z);
 free(ifixed);
 free(ifree);
 free(dt);
 free(b);
 free(temptau);
 free(ik);
}
