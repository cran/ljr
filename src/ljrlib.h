#include <stdlib.h>
#include <time.h>
#include <math.h>
#define repeat for(;;)

// Calculates x=x+dx.
// x and dx are n-dimensional vectors.
void radd(double *x,double *dx,int n);

// Calculates x=x-dx.
// x and dx are n-dimensional vectors.
void rsub(double *x,double *dx,int n);

// Calculates the sup norm for a vector.
// x is an n-dimensional vectors.
// res is the result
void rdif(double *x,int N,double *res);

// Computes inverse logit transformation.
// eta is a vector of N logits.
// The resulting vector of inverses is p=e^eta/(1+e^eta).
void ril(double *eta,double *p,int N);

// This function computes the matrix product A*W*t(A).
// A is a p-by-n matrix by row.
// w is a vector containing the elements on the diagonal of
// a n-by-n diagonal matrix W.
// The result is stored in R with the entries stored by rows.
void rmdm(double *A,double *w,double *R,int p,int n);

// This function multiplies a matrix by a vector.
// A is a m-by-n matrix; its entries are stored by rows.
// b is a n-dimensional vector.
// The result is stored in r with entries stored by rows.
void rmv(double *A, double *b, double *r, int m, int n);

// This function transposes a matrix.
// A is a m-by-n matrix by row.
// The resulting transpose is the n-by-m matrix B.
void rt(double *A,double *B, int m, int n);

// Computes weights based on probabilities and sizes.
// p is the vector of probabilities.
// n is the vector of sizes.
// N is the dimension of all vectors.
// w is the resulting vector of weights n*p*(1-p).
void rw(double *p,double *n,double *w,int N);

// Computes y-n*p.
// y is the vector of binomial responses.
// p is the vector of probabilities.
// n is the vector of sizes.
// N is the dimension of all vectors.
void rynp(double *y,double *n,double *p,double *res,int N);

// Computes empirical estimates of p.
// y is the vector of binomial responses.
// p is the vector of probabilities.
// n is the vector of sizes.
// N is the dimension of all vectors.
void rem(double *y,double *n,double *p,int N);

// Computes transformed y used for initial fitting to
// obtain starting point for Newton-Raphson.
// Specifically, newy=log(p/(1-p)).
// p is the vector of sizes.
// newy is the transformed y used to obtain starting point for NR.
void rny(double *p,double *newy,int N);

// Computes the rescaled log-likelihood function for fitted etas.
// y is the vector of binomial responses.
// n is the vector of sizes.
// eta is the fitted logit at the fitted responses.
// res is the value of the log-likelihood function corresponding to the estimates
// N is the dimension of all vectors.
void rlk(double *y,double *n,double *eta,double *res,int N);

// Fits logistic regression with starting point based on linear regression.
// y is the N-dimensional vector of binomial responses.
// n is the N-dimensional vector of sizes.
// x is the m-by-N matrix of explanatory variables, stored by row.
// ofst is the N-dimensional vector of offsets.
// b is the m-dimensional vector of coefficient estimates computed
// by the algorithm.
// wlik is the scalar fitted log-likelihood (stored as pointer for use in R).
void lr(double *y,double *n,double *x,double *ofst,double *b,int N,int m,double *wlik);

// Builds design matrix for joinpoint model.
// z is the design matrix for the joinpoint algorithm.
// x is the m-by-N matrix of explanatory variables, stored by row.
// tm is the N-dimensional vector of observation times.
// ifixed is a vector of the indices of the fixed joinpoints.
// ifree is a vector of the indices of the free joinpoints.
// ncps is a 2-dimensional vector with the number of fixed and free 
// joinpoints, respectively.
// newz>0 indicates that z already contains x in the first m rows.
void rz(double *z,double *x,double *tm,int *ifixed,int *ifree,int *ncps,int N,int m,int newz);

// Computes joinpoint model estimates based on lr fit.
// beta is a vector of intercept and slope coefficients
// gamma is a vector of coefficients corresponding to tm and the terms with 
// joinpoints.
// b is the (m+1)-dimensional vector of coefficient estimates computed from lr.
// ncps is a 2-dimensional vector with the number of fixed and free 
// joinpoints, respectively.
void rcf(double *beta,double *gamma,double *b,int m,int *ncps);

// Fits logistic joinpoint regression with no joinpoints.
// y is the N-dimensional vector of binomial responses.
// n is the N-dimensional vector of sizes.
// tm is the N-dimensional vector of observation times.
// x is the m-by-N matrix of explanatory variables, stored by row.
// ofst is the N-dimensional vector of offsets.
// beta is the (m+1)-dimensional vector of coefficient estimates computed 
// by the algorithm.
// gamma is the estimate of the coefficient for the term with the times.
// zlik is the scalar fitted log-likelihood.
void ljr0(double *y,double *n,double *tm,double *x,double *ofst,double *beta,double *gamma,int *Nptr,int *mptr,double *zlik);

// Fits logistic joinpoint regression with 1 joinpoint.
// y is the N-dimensional vector of binomial responses.
// n is the N-dimensional vector of sizes.
// tm is the N-dimensional vector of observation times.
// x is the m-by-N matrix of explanatory variables, stored by row.
// ofst is the N-dimensional vector of offsets.
// beta is the (m+1)-dimensional vector of coefficient estimates computed 
// by the algorithm.
// gamma is the estimate of the coefficient for the term with the times.
// tau is the estimate of the joinpoint.
// zlik is the scalar fitted log-likelihood.
void ljr1(double *y,double *n,double *tm,double *x,double *ofst,double *beta,double *gamma,double *tau,int *Nptr,int *mptr,double *zlik);

// Fits logistic joinpoint regression with 2 joinpoints.
// y is the N-dimensional vector of binomial responses.
// n is the N-dimensional vector of sizes.
// tm is the N-dimensional vector of observation times.
// x is the m-by-N matrix of explanatory variables, stored by row.
// ofst is the N-dimensional vector of offsets.
// beta is the (m+1)-dimensional vector of coefficient estimates computed 
// by the algorithm.
// gamma is the estimate of the coefficient for the term with the times.
// tau contains the estimates of the joinpoints.
// zlik is the scalar fitted log-likelihood.
void ljr2(double *y,double *n,double *tm,double *x,double *ofst,double *beta,double *gamma,double *tau,int *Nptr,int *mptr,double *zlik);

double fmin2(double x, double y);

double fmax2(double x, double y);

void rrbinom(int *nin, double *pin, double *x);

// This function generates binomial random variables according to the logistic
// joinpoint regression model.
// beta is the (m+1)-dimensional vector of coefficients corresponding to the
// variables in x.
// gamma is the vector of coefficient estimates corresponding to tm and the
// corresponding change points.
// tau is the vector of joinpoints.
// n is the N-dimensional vector of sizes.
// tm is the N-dimensional vector of observation times.
// x is the m-by-N matrix of explanatory variables, stored by row.
// ofst is the N-dimensional vector of offsets.
// y is the N-dimensional vector of binomial responses.
void rgy(double *beta,double *gamma,double *tau,double *n,double *tm,double *x,double *ofst,double *y,int N,int m,int ncps);

// This function tests the null hypothesis of 0 joinpoints vs the alternative
// that there is joinpoint via Monte Carlo simulations.
// y is the N-dimensional vector of binomial responses.
// n is the N-dimensional vector of sizes.
// tm is the N-dimensional vector of observation times.
// x is the m-by-N matrix of explanatory variables, stored by row.
// ofst is the N-dimensional vector of offsets.
// R is the number of Monte Carlo simulations.
// p is the p-value returned by the function.
void test01(double *y,double *n,double *tm,double *x,double *ofst,int *N,int *m,int *R,double *p);

// This function performs the backward algorithm for logistic regression.
// y is the N-dimensional vector of binomial responses.
// n is the N-dimensional vector of sizes.
// tm is the N-dimensional vector of observation times.
// x is the m-by-N matrix of explanatory variables, stored by row.
// ofst is the N-dimensional vector of offsets.
// R is the number of Monte Carlo simulations for each test.
// p is the vector of p-values for the tests.
// alpha is the size of each test.
// nulls is the number of joinpoints under the null hypotheses.
// alts is the number of joinpoints under the alt hypotheses.
void backward2(double *y,double *n,double *tm,double *x,double *ofst,int *N,int *m,int *R,double *p,double *alpha,int *nulls,int *alts);

// This function performs the forward algorithm for logistic regression.
// y is the N-dimensional vector of binomial responses.
// n is the N-dimensional vector of sizes.
// tm is the N-dimensional vector of observation times.
// x is the m-by-N matrix of explanatory variables, stored by row.
// ofst is the N-dimensional vector of offsets.
// R is the number of Monte Carlo simulations for each test.
// p is the vector of p-values for the tests.
// alpha is the size of each test.
// nulls is the number of joinpoints under the null hypotheses.
// alts is the number of joinpoints under the alt hypotheses.
void forward2(double *y,double *n,double *tm,double *x,double *ofst,int *N,int *m,int *R,double *p,double *alpha,int *nulls,int *alts);

void rgi(double *tm,int *ti,int *numti,int n);

void mvfixed(double *z,double *tm,int *ncps,int *dt,int n,int m,int v,int loc);

void mvfree(double *z,double *tm,int *dt,int n,int m,int v,int loc);

void fillfixed(double *z,double *tm,int *ncps,int *dt,int n,int m,int v,int loc,int skp,int upto);

void fillfree(double *z,double *tm,int *ncps,int *dt,int n,int m,int *ifree);

// Removes the kth row of an p-by-N matrix z.
void rzrmrow(double *z,int k,int N,int p);

void mvfixedrm(double *z,double *tm,int *ncps,int *dt,int n,int m,int v,int loc);

void mvfreermint(double *z,double *tm,int *dt,int n,int m,int v,int loc);

void mvfreermtm(double *z,double *tm,int *dt,int n,int m,int v,int loc);

void rcfrmint(double *beta,double *gamma,double *b,int m,int *ncps);

void rcfrmtm(double *beta,double *gamma,double *b,int m,int *ncps);

void ljr1rmint(double *y,double *n,double *tm,double *x,double *ofst,double *beta,double *gamma,double *tau,int *Nptr,int *mptr,double *zlik);

void ljr1rmtm(double *y,double *n,double *tm,double *x,double *ofst,double *beta,double *gamma,double *tau,int *Nptr,int *mptr,double *zlik);

void ljr11(double *y,double *n,double *tm,double *x,double *ofst,int *Nptr,int *mptr,int *Rptr,double *p);

void fillfixedrm(double *z,double *tm,int *ncps,int *dt,int n,int m,int v,int loc,int skp
,int upto);

void fillfreermint(double *z,double *tm,int *ncps,int *dt,int n,int m,int *ifree);

void fillfreermtm(double *z,double *tm,int *ncps,int *dt,int n,int m,int *ifree);

void ljr2rmint(double *y,double *n,double *tm,double *x,double *ofst,double *beta,double *gamma,double *tau,int *Nptr,int *mptr,double *zlik);

void ljr2rmtm(double *y,double *n,double *tm,double *x,double *ofst,double *beta,double *gamma,double *tau,int *Nptr,int *mptr,double *zlik);

void ljr22(double *y,double *n,double *tm,double *x,double *ofst,int *Nptr,int *mptr,int *Rptr,double *p);

void test02(double *y,double *n,double *tm,double *x,double *ofst,int *N,int *m,int *R,double *p);

void test12(double *y,double *n,double *tm,double *x,double *ofst,int *N,int *m,int *R,double *p);

int checkseq(int *x,int M,int k,int *numfixed);

void srz(int *ik,int *dt,int *ifixed,int *ifree,int *ncps,int k,int M,int nfixed);

void ljrk(int *kptr,double *y,double *n,double *tm,double *x,double *ofst,double *beta,double *gamma,double *tau,int *Nptr,int *mptr,double *zlik);

void backwardk(int *Kptr,double *y,double *n,double *tm,double *x,double *ofst,int *Nptr,int *mptr,int *Rptr,double *p,double *alphaptr,int *nulls,int *alts);

void forwardk(double *y,double *n,double *tm,double *x,double *ofst,int *Nptr,int *mptr,int *Rptr,double *p,double *alphaptr,int *ncps);

void testjk(int *jptr,int *kptr,double *y,double *n,double *tm,double *x,double *ofst,int *Nptr,int *mptr,int *Rptr,double *p);

void ljrkrmint(int *kptr,double *y,double *n,double *tm,double *x,double *ofst,double *beta,double *gamma,double *tau,int *Nptr,int *mptr,double *zlik);

void ljrkrmtm(int *kptr,double *y,double *n,double *tm,double *x,double *ofst,double *beta,double *gamma,double *tau,int *Nptr,int *mptr,double *zlik);

void ljrkk(int *kptr,double *y,double *n,double *tm,double *x,double *ofst,int *Nptr,int *mptr,int *Rptr,double *p);
