void rcf(double *beta,double *gamma,double *b,int m,int *ncps){
 int i;
 for (i=1;i<m+1;i++)
  beta[i]=b[i-1];
 beta[0]=b[m];
 gamma[0]=b[m+1];
 if (ncps[1]>0)
  for (i=1;i<=ncps[1];i++)
   gamma[i]=b[m+2*i+1]-b[m+2*i-1];
 if (ncps[0]>0)
  for (i=1;i<=ncps[0];i++)
   gamma[ncps[1]+i]=b[m+2*ncps[1]+1+i];
}
