void rw(double *p,double *n,double *w,int N){
 int i;
 for (i=0;i<N;i++)
  w[i]=n[i]*p[i]*(1-p[i]);
}
