void rem(double *y,double *n,double *p,int N){
 int i;
 for (i=0;i<N;i++){
  p[i]=y[i]/n[i];
  if (p[i]==0)
   p[i]=.5/(n[i]+.5);
  if (p[i]==1)
   p[i]=n[i]/(n[i]+.5);
 }
}
