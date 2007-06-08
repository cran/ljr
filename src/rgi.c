void rgi(double *tm,int *ti,int *numti,int n){
 int i,count=1;
 ti[0]=0;
 for (i=1;i<n;i++)
  if (tm[i]>tm[i-1])
   ti[count++]=i;
 ti[count]=n;
 *numti=count;
}
