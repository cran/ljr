#include "ljrlib.h"
int checkseq(int *x,int M,int k,int *numfixed){
 int fixed=0;
 int i,j;
 int res=1;
 while ((fixed<k)&&(x[fixed]<=M)){
  fixed++;
 }
 for (j=fixed;j<k;j++){
  for (i=0;i<fixed;i++)
   if ((x[i]==(x[j]%M))|(x[i]==(x[j]%M)+1))
    goto finis1;
   if (j<k-1)
    if (x[j]==x[j+1]-1)
     goto finis1;
   if (j>0)
    if ((x[j]==x[j-1]+1))
     goto finis1;
 }
 if (res==1)
  goto finis2;
 finis1:
  res=0;
  goto finis2;
 finis2:
 *numfixed=fixed;
 return res;
}
