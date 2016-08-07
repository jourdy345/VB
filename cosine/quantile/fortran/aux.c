#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

void F77_SUB(rndstart)(void)
{
  GetRNGstate();
}

void F77_SUB(rndend)(void)
{
  PutRNGstate();
}

double F77_SUB(rndunif)(void)
{
  return unif_rand();
}

double F77_SUB(rndnorm)(void)
{
  return norm_rand();
}

double F77_SUB(rvunif)(double *a, double *b)
{
  return runif(*a, *b);
}

double F77_SUB(rvnorm)(double *mu, double *sigma)
{
  return rnorm(*mu, *sigma);
}

double F77_SUB(rvgamma)(double *shape, double *scale)
{
  return rgamma(*shape, *scale);
}

double F77_SUB(rvbeta)(double *a, double *b)
{
  return rbeta(*a, *b);
}

double F77_SUB(rvexp)(double *scale)
{
  return rexp(*scale);
}