#include <R.h>

void prpm(char name, double *A, int m, int n){
 int i,j;
 Rprintf("%c=\n",name);
 for (i=0;i<m;i++){
  for (j=0;j<n;j++)
   Rprintf("%10.4f ",*(A++));
  Rprintf("\n");
 }
}
