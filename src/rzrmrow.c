void rzrmrow(double *z,int k,int N,int p){
 int a,i; 
 a=N*(k-1);
 for (i=0;i<N*(p-k);i++){
  z[a]=z[a+N];
  a++;
 }
}
