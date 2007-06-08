void fillfixed(double *z,double *tm,int *ncps,int *dt,int n,int m,int v,int loc,int skp,int upto){
 int i,a,dto,dtn,dte;
 dto=dt[loc];
 dtn=dt[skp];
 dte=dt[upto];
 if (dtn<dto){
  a=n*(m+v+2*(ncps[1]==0)+2*ncps[1]-1)+dtn;
  for (i=dtn;i<dto;i++)
   z[a++]=0;
  for (i=dto;i<dte;i++)
   z[a++]=tm[i]-tm[dto-1];
 }
 else{
  a=n*(m+v+2*(ncps[1]==0)+2*ncps[1]-1)+dtn;
  for (i=dtn;i<dte;i++)
   if (dto>0)
    z[a++]=tm[i]-tm[dto-1];
   else
    z[a++]=tm[i];
 }
}
