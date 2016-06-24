// [[Rcpp::depends(RcppArmadillo)]]
#define _USE_MATH_DEFINES
#include <RcppArmadillo.h>
#include <cmath>
#include <Rmath.h>
#include <random>
using namespace Rcpp; using namespace arma;


// [[Rcpp::export]]
arma::cube setVarPhi(arma::vec x, int J) {
  int n = x.size();
  arma::cube varphi(J+1, J+1, n);
  for (int i = 0; i < n; ++i) {
    arma::mat temp(J+1, J+1);
    temp(0,0) = x(i) - 0.5;
    
    for (int j = 1; j < J+1; ++j) {
      
      temp(0, j) = std::sqrt(2.0)/(M_PI * double(j)) * std::sin(M_PI * double(j) * x(i)) - std::sqrt(2.0) / std::pow((M_PI * double(j)), 2.0) * (1.0 - std::cos(M_PI * double(j)));
      temp(j, 0) = temp(0, j);
      
      for (int k = j; k < J+1; ++k) {
        
        if (k == j) {
          temp(j, j) = std::sin(2.0 * M_PI * double(j) * x(i)) / (2.0 * M_PI * double(j)) + x(i) - 0.5;
        } else {
          temp(j, k) = std::sin(M_PI * (double(j + k)) * x(i)) / (M_PI * (double(j + k))) + std::sin(M_PI * (double(j - k)) * x(i)) / (M_PI * (double(j - k))) - (1.0 - std::cos(M_PI * (double(j + k)))) / std::pow((M_PI * (double(j + k))), 2.0) - (1.0 - std::cos(M_PI * (double(j - k)))) / std::pow((M_PI * (double(j - k))), 2.0);
          temp(k, j) = temp(j, k);
        }
      }
    }
    varphi.slice(i) = temp;
  }
  return varphi;
}