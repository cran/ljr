#include <R.h>

void iprpm(char name, int *A, int m, int n){
 int i,j;
 Rprintf("%c=\n",name);
 for (i=0;i<m;i++){
  for (j=0;j<n;j++)
   Rprintf("%i ",*(A++));
  Rprintf("\n");
 }
}
