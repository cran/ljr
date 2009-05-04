#include <math.h>
void rny(double *p,double *newy,int N){
 int i;
 for (i=0;i<N;i++)
  newy[i]=log(p[i]/(1-p[i]));
}
