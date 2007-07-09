#include <R.h>
#include <Rmath.h>

void rrbinom(int *nin, double *pin, double *x)
{
    int n = nin[0];
    double p = pin[0];
    GetRNGstate();
    x[0]=0.0+rbinom(n, p);
    PutRNGstate();
}
