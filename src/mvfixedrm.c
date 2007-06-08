void mvfixedrm(double *z,double *tm,int *ncps,int *dt,int n,int m,int v,int loc){
 int i,a,dto,dtn;
 dto=dt[loc-1];
 dtn=dt[loc];
 a=n*(m+v+2*ncps[1])+dto;
 for (i=dto;i<dtn;i++)
  z[a++]=0;
 for (i=dtn;i<n;i++)
  z[a++]=tm[i]-tm[dto];
}
