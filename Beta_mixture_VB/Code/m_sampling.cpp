// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <Rmath.h>
using namespace Rcpp; using namespace arma;

// [[Rcpp:export]]
double my_fun(double x, double Cm, double a, double b) {
  return exp(-Cm*x) * pow(x, a-1) * pow(1-x, b-1)
}

// [[Rcpp::export]]
double opt(Function f) {

}