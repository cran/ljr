void rmdm(double *A,double *w,double *R,int p,int n){
 int i=0,j,k,m=0;
 int h1,h2=0;
 int d,r;
/* Initialize results to 0. */
 for (j=0;j<p*p;j++)
  R[j]=0;
/* Start filling in first row. */
 for (j=0;j<p;j++){
  h1=0;
  for (k=0;k<n;k++)
   R[m]+=A[h1++]*A[h2++]*w[k];
  m++;
 }
 for (i=1;i<p;i++){
/* Use symmetry for lower triangle. */
  for (j=0;j<i;j++){
   d=m/p;
   r=m-d*p;
   R[m]=R[p*r+d];
   m++;
  }
/* Fill in the rest of the row. */
  h2=i*n;
  for (j=i;j<p;j++){
   h1=i*n;
   for (k=0;k<n;k++)
    R[m]+=A[h1++]*A[h2++]*w[k];
   m++;
  }
 }
}
