void rz(double *z,double *x,double *tm,int *ifixed,int *ifree,int *ncps,int N,
		int m,int newz){
 int i,j,k;
 int a;
 if (newz){
  a=0;
  for (i=0;i<(m*N);i++)
   z[a++]=x[i];
 }
 else 
  a=m*N;
 if (ncps[1]==0){
  for (i=0;i<N;i++)
   z[a++]=1;
  for (i=0;i<N;i++)
   z[a++]=tm[i];
 }
 else{
  for (j=0;j<ifree[0];j++)
   z[a++]=1;
  for (j=ifree[0];j<N;j++)
   z[a++]=0;
  for (j=0;j<ifree[0];j++)
   z[a++]=tm[j];
  for (j=ifree[0];j<N;j++)
   z[a++]=0;
  for (k=1;k<ncps[1];k++){
   for (j=0;j<ifree[k-1];j++)
    z[a++]=0;
   for (j=ifree[k-1];j<ifree[k];j++)
    z[a++]=1;
   for (j=ifree[k];j<N;j++)
    z[a++]=0;
   for (j=0;j<ifree[k-1];j++)
    z[a++]=0;
   for (j=ifree[k-1];j<ifree[k];j++)
    z[a++]=tm[j];
   for (j=ifree[k];j<N;j++)
    z[a++]=0;
  }
  for (j=0;j<ifree[ncps[1]-1];j++)
   z[a++]=0;
  for (j=ifree[ncps[1]-1];j<N;j++)
   z[a++]=1;
  for (j=0;j<ifree[ncps[1]-1];j++)
   z[a++]=0;
  for (j=ifree[ncps[1]-1];j<N;j++)
   z[a++]=tm[j];
 }
 for (k=0;k<ncps[0];k++){
  for (j=0;j<ifixed[k];j++)
   z[a++]=0;
  for (j=ifixed[k];j<N;j++)
   z[a++]=tm[j]-tm[ifixed[k]-1];
 }
}
