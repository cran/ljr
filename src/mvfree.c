void mvfree(double *z,double *tm,int *dt,int n,int m,int v,int loc){
 int i,ia,a,dto,dtn;
 dto=dt[loc-1];
 dtn=dt[loc];
 ia=n*(m+2*v-2);
 a=ia+dto;
 for (i=dto;i<dtn;i++)
  z[a++]=1;
 a=ia+n+dto;
 for (i=dto;i<dtn;i++)
  z[a++]=tm[i];
  a=ia+2*n+dto;
 for (i=dto;i<dtn;i++)
  z[a++]=0;
 a=ia+3*n+dto;
 for (i=dto;i<dtn;i++)
  z[a++]=0;
}
