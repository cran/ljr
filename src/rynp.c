void rynp(double *y,double *n,double *p,double *res,int N){
 int i;
 for (i=0;i<N;i++)
  res[i]=y[i]-n[i]*p[i];
}
