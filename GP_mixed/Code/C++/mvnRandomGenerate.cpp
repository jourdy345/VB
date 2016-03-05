#include <iostream>
#include <cmath>
#include <armadillo>
#include "mvnRandomGenerate.h"


arma::mat mvnRandomGenerate(int n, arma::vec mu, arma::mat sigma) {
   int nCols   = sigma.n_cols;
   arma::mat Y = arma::randn(n, nCols);

   return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}