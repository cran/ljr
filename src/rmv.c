void rmv(double *A, double *b, double *r, int m, int n){
  int i,j;
  double v,*pb;
  for (i=0;i<m;i++){
   v=0;
   pb=b;
   for (j=0;j<n;j++)
    v+=(*(A++))*(*(pb++));
   *(r++)=v;
  }
}
