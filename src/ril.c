#include <math.h>
void ril(double *eta,double *p,int N){
 int i;
 for (i=0;i<N;i++)
  p[i]=1-1/(1+exp(eta[i]));
}
