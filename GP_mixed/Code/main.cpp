#define _USE_MATH_DEFINES
#include <iostream>
#include <armadillo>
#include <cmath> //for erfc function and math constants

int main() {
  arma::mat x(6, 6); x.randn();
  arma::mat y(6, 6); y.randn();
  arma::rowvec z = x.row(2) % y.row(3);
  arma::colvec zt = z.t();
  std::cout << x << std::endl;
  std::cout << y << std::endl;
  std::cout << zt << std::endl;
  return 0;
}