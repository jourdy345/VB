// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <Rmath.h>
using namespace Rcpp; using namespace arma;

// [[Rcpp::export]]
arma::mat blockdiag( arma::field<mat> x ) {

  //x: list of matrices 

  unsigned int n = x.n_rows;
  int dimen = 0;
  arma::ivec dimvec(n);

  for(unsigned int i=0; i<n; i++) {
      dimvec(i) = x(i,0).n_rows; 
      dimen += dimvec(i);
  }

  mat X(dimen,dimen,fill::zeros);
  int idx=0;

  for(unsigned int i=0; i<n; i++) {
      X.submat( idx, idx, idx + dimvec(i) - 1, idx + dimvec(i) - 1 ) = x(i,0);
      idx = idx + dimvec(i);
  }

  return(X);
}