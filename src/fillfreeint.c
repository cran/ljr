void fillfreermint(double *z,double *tm,int *ncps,int *dt,int n,int m,int *ifree){
 int j,k,a;
 a=m*n;
 for (j=0;j<ifree[0];j++)
  z[a++]=tm[j];
 for (j=ifree[0];j<n;j++)
  z[a++]=0;
 for (k=1;k<ncps[1];k++){
  for (j=0;j<ifree[k-1];j++)
   z[a++]=0;
  for (j=ifree[k-1];j<ifree[k];j++)
   z[a++]=1;
  for (j=ifree[k];j<n;j++)
   z[a++]=0;
  for (j=0;j<ifree[k-1];j++)
   z[a++]=0;
  for (j=ifree[k-1];j<ifree[k];j++)
   z[a++]=tm[j];
  for (j=ifree[k];j<n;j++)
   z[a++]=0;
 }
 for (j=0;j<ifree[ncps[1]-1];j++)
  z[a++]=0;
 for (j=ifree[ncps[1]-1];j<n;j++)
  z[a++]=1;
 for (j=0;j<ifree[ncps[1]-1];j++)
  z[a++]=0;
 for (j=ifree[ncps[1]-1];j<n;j++)
  z[a++]=tm[j];
}
