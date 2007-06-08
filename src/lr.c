#include "ljrlib.h"
void lr(double *y,double *n,double *x,double *ofst,double *b,int N,int m,double *wlik){
  char uplo = 'U';
  int o=1;
  int info;
  int maxit=20;
  int count=1;
  int i,j;

  double *db;
  double *eta;
  double *w;
  double *p;
  double *yerr;
  double *tx;
  double *xwx;

  double ep=1e-7;
  double diff[1];

  db=malloc(m*sizeof(double));
  eta=malloc(N*sizeof(double));
  w=malloc(N*sizeof(double));
  p=malloc(N*sizeof(double));
  yerr=malloc(N*sizeof(double));
  tx=malloc(N*m*sizeof(double));
  xwx=malloc(m*m*sizeof(double));

// Pick starting point by performing linear regression on the logit.
  rt(x,tx,N,m);
  rem(y,n,p,N);
  rny(p,eta,N);
  rsub(eta,ofst,N);
  rw(p,n,w,N);
  rmdm(x,w,xwx,m,N);
  for (i=0;i<N;i++)
   yerr[i]=w[i]*eta[i];
  rmv(x,yerr,b,m,N);
  dposv_(&uplo,&m,&o,xwx,&m,b,&m,&info);
  *diff=ep+1;

// Iterative re-weighted least squares procedure.
  while ((*diff>ep)&(count<maxit)){
   rmv(tx,b,eta,N,m);
   radd(eta,ofst,N);
   ril(eta,p,N);
   rw(p,n,w,N);
   rmdm(x,w,xwx,m,N);
   rynp(y,n,p,yerr,N);
   rmv(x,yerr,db,m,N);
   dposv_(&uplo, &m, &o, xwx, &m, db, &m, &info);
   rdif(db,m,diff);
   radd(b,db,m);     
   count++;
  }

// If the above algorithm does not converge, revert to 
// an iterative grid search to choose the starting point.
  if (*diff>ep){
   double oldcoef;
   double tempb[m];
   double tlik[1];
   *tlik=0;
   for (i=0;i<m;i++){
    b[i]=0;
    tempb[i]=0;
   }
   rt(x,tx,N,m); 
   rmv(tx,b,eta,N,m);
   radd(eta,ofst,N);
   rlk(y,n,eta,wlik,N);
   *diff=1;
   while (*diff==1){
    *diff=0;  
    for (i=0;i<m;i++){
     oldcoef=b[i];
     for (j=-100;j<=100;j++){
      tempb[i]=j/100.;
      rmv(tx,tempb,eta,N,m);
      radd(eta,ofst,N);
      rlk(y,n,eta,tlik,N);
      if (*tlik>*wlik){
       *wlik=*tlik;
       b[i]=j/100.;  
      }
      if ((b[i]==1)||(b[i]==-1)){
       for (j=-100;j<=100;j++){
        tempb[i]=j/10.;
        rmv(tx,tempb,eta,N,m);
        radd(eta,ofst,N);
        rlk(y,n,eta,tlik,N);
        if (*tlik>*wlik){
         *wlik=*tlik;
         b[i]=j/10.;  
        }   
        if ((b[i]==.1)||(b[i]==-.1)){
         for (j=-100;j<=100;j++){
          tempb[i]=j;
          rmv(tx,tempb,eta,N,m);
          radd(eta,ofst,N);
          rlk(y,n,eta,tlik,N);
          if (*tlik>*wlik){
           *wlik=*tlik;
           b[i]=j;  
          }
         }
        }
       }
      }
     }
     if (b[i]!=oldcoef)
      *diff=1;
     tempb[i]=b[i];
    }
   }

// Iterative re-weighted least squares procedure.
   *diff=ep+1;   
   count=1;
   while ((*diff>ep)&(count<maxit)){
    rmv(tx,b,eta,N,m);
    radd(eta,ofst,N);
    ril(eta,p,N);
    rw(p,n,w,N);
    rmdm(x,w,xwx,m,N);
    rynp(y,n,p,yerr,N);
    rmv(x,yerr,db,m,N);
    dposv_(&uplo, &m, &o, xwx, &m, db, &m, &info);
    rdif(db,m,diff);
    radd(b,db,m);     
    count++;
   }
  }
 rlk(y,n,eta,wlik,N);
 free(db);
 free(eta);
 free(w);
 free(p);
 free(yerr);
 free(tx);
 free(xwx);
}
