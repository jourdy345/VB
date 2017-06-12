// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <Rmath.h>
using namespace Rcpp; using namespace arma;

// [[Rcpp::export]]
arma::mat duplicate_matrix (int n) {
  arma::mat mat1 = arma::eye<mat>(n, n);
  arma::vec index = arma::vectorise(arma::linspace<vec>(1, n*(n+1)/2, n*(n+1)/2));
  for (int j = 0; j < n; ++j) {
    for (int i = j; i < n; ++i) {
      if (j == 0) mat1(i, j) = index(i);
      else {
        mat1(i, j) = index(n*j - (j-1)*j/2 + i - j, 0);
      }
    }
  }
  arma::mat mat2 = arma::symmatl(mat1);
  arma::vec temp_vec = arma::vectorise(mat2);
  int t = temp_vec.size();
  int s = index.size();
  arma::mat result(t, s);
  for (int k = 0; k < t; ++k) {
    for (int u = 0; u < s; ++u) {
      if (temp_vec(k) == index(u)) result(k, u) = 1;
      else result(k, u) = 0;
    }
  }
  std::cout << mat1 << std::endl;
  return result;
}