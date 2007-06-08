#include <math.h>
void rlk(double *y,double *n,double *eta,double *res,int N){
 int i;
 double mn;
 mn=n[0];
 for (i=1;i<N;i++)
  if (n[i]>mn)
   mn=n[i];
 *res=0;
 for (i=0;i<N;i++)
  *res+=(eta[i]*y[i]-n[i]*log(1+exp(eta[i])))/mn; 
}
