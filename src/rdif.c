#include <math.h>
void rdif(double *x,int N,double *res){
 double d;
 int i;
 *res=0;
 for (i=0;i<N;i++){ 
  d=fabs(x[i]);
  if (d>*res)
   *res=d;
 }
}
