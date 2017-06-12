#include "./math/stan/math.hpp"
#include <iostream>
#include <vector>
#include <random>
using Eigen::Matrix;
using Eigen::Dynamic;

// struct normal_ll {
//   const Matrix<double,Dynamic,1> y_;

//   normal_ll(const Matrix<double,Dynamic,1>& y) : y_(y) {}
//   template <typename T>
//   T operator() (const Matrix<T,Dynamic,1>& theta) const {
//     T mu = theta[0];
//     T sigma = theta[1];
//     T lp = 0;
//     for (int n = 0; n < y_.size(); ++n) {
//       lp += stan::math::normal_log(y_[n],mu,sigma);
//     }
//     return lp;
//   }
// };

int main() {
  
  // set RNG
  std::random_device rd;
  std::mt19937 gen(rd());

  std::normal_distribution<double> randn(0.,1.);

  Matrix<double,Dynamic,1> X(3);
  Matrix<double,Dynamic,1> X_(3);
  for (int i = 0; i < 3; ++i) {
    X(i) = randn(gen);
    X_(i) = randn(gen);
  }
  Matrix<double,1,Dynamic> Xt = X_.transpose();
  std::cout << stan::math::multiply(X,Xt) << std::endl;
  std::cout << X*X_.transpose() << std::endl;
  // return 0;
}