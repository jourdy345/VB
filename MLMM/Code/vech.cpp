// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <Rmath.h>
using namespace Rcpp; using namespace arma;

// [[Rcpp::export]]
arma::vec vech(arma::mat X) {
  int n = X.n_rows;
  arma::vec temp(n*(n+1)/2);
  for (int j = 0; j < n; ++j) {
    for (int i = j; i < n; ++i) {
      temp(n*j - (j-1)*j/2 + i - j) = X(i, j);
    }
  }
  return temp;
}