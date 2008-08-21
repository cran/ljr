void rt(double *A,double *B, int m, int n){
 int i,j;
 int p=m*n-1;
 for (i=0;i<n;i++){
  for (j=0;j<m;j++){
   *B=*(A++);  
   B+=n;
  }
  B-=p;
 }
}
