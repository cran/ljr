#include "ljrlib.h"
void srz(int *ik,int *dt,int *ifixed,int *ifree,int *ncps,int k,int M,int nfixed){
 int i;
 ncps[0]=nfixed;
 ncps[1]=k-nfixed;
 for (i=0;i<nfixed;i++)
  ifixed[i]=dt[ik[i]+1];
 for (i=0;i<k-nfixed;i++)
  ifree[i]=dt[ik[i+nfixed]%M+1];
}
