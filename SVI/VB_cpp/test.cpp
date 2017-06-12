#include "./math/stan/math.hpp"
#include <iostream>
#include <vector>
using Eigen::Matrix;
using Eigen::Dynamic;
// #include "logh.hpp"
/* sudo g++ -I ~/Desktop/tmp/VB_cpp/math -I ~/Desktop/tmp/VB_cpp/math/lib/eigen_3.2.9/ -I ~/Desktop/tmp/VB_cpp/math/lib/boost_1.62.0/ -I ~/Desktop/tmp/VB_cpp/math/lib/cvodes_2.9.0/include/ -std=c++11 -stdlib=libc++ test.cpp -o a.out */
/* sudo clang++ -std=gnu++11 -stdlib=libc++ -I ~/Desktop/tmp/VB_cpp/math -I ~/Desktop/tmp/VB_cpp/math/lib/eigen_3.2.9/ -I ~/Desktop/tmp/VB_cpp/math/lib/boost_1.62.0/ -I ~/Desktop/tmp/VB_cpp/math/lib/cvodes_2.9.0/include/ test.cpp -o a.out */




// template <typename T>
// T test_fun(const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>& X, const Eigen::Matrix<T,Eigen::Dynamic,1>& beta) {
//   double res = stan::math::dot_product(beta,(X*beta));
//   return res;
// }












// template <typename T>
// T logh(const Eigen::Matrix<T,Eigen::Dynamic,1> parameters, const Eigen::Matrix<T,Eigen::Dynamic,1> hyperparameters,
//        const Eigen::Matrix<T,Eigen::Dynamic,1> y, const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> vphi,
//        const int n, const int J, const int B)
// {
//   Eigen::MatrixXd theta(J,1);
//   for (int i = 0; i < J; ++i) {
//     theta(i) = parameters(i);
//   }
//   T psi = parameters(J);
//   T sig2 = parameters(J+1);
//   T tau2 = parameters(J+2);
//   T gam  = parameters(J+3);

//   // extract hyperparameters
//   T r_sig = hyperparameters(0);
//   T s_sig = hyperparameters(1);
//   T r_tau = hyperparameters(2);
//   T s_tau = hyperparameters(3);
//   T a     = hyperparameters(4);
//   T b     = hyperparameters(5);
//   T w0    = hyperparameters(6);
//   Eigen::MatrixXd diagm(B,B);
//   for (int i=0; i<B;++i) {
//     for (int j=0;j<B;++j) {
//       if (i==j) {
//         diagm(i,i) = 1.;
//       } else {
//         diagm(i,j) = 0.;
//       }
//     }
//   }
//   Eigen::VectorXd mu         = vphi*theta;
//   Eigen::MatrixXd Covariance = (sig2*diagm.array()).matrix();
//   double res = stan::math::multi_normal_log(y,mu,Covariance)/double(B)*double(n)+stan::math::inv_gamma_log(sig2,r_sig/2.,s_sig/2.)+stan::math::inv_gamma_log(tau2,r_tau/2.,s_tau/2.)+stan::math::inv_gamma_log(psi,a,b)+stan::math::log(w0)-w0*gam;
//   for (int j = 0; j < J; ++j) {
//     res += stan::math::normal_log(theta(j),0.,sig2*tau2*stan::math::exp(-(double(j)+1.)*gam));
//   }
//   return res;
// }

// template logh<double>(const Eigen::Matrix<double,Eigen::Dynamic,1> parameters, const Eigen::Matrix<double,Eigen::Dynamic,1> hyperparameters, const Eigen::Matrix<double,Eigen::Dynamic,1> y, const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> vphi, const int n, const int J, const int B);
struct quadratic {
  // const Eigen::Matrix<double,Eigen::Dynamic,1> beta_;
  // const Eigen::MatrixXd X_ = X;
  const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> X_;

  quadratic(const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>& X) : X_(X) {}
  template <typename T>
  T operator()(const Eigen::Matrix<T,Eigen::Dynamic,1>& beta) const {
    T res = 0;
    Eigen::Matrix<T,Dynamic,1> tmp;
    tmp = stan::math::multiply(X_,beta);
    res = stan::math::dot_product(beta,tmp);
    return res;
  }
};

int main() {
  // stan::math::var res = 0.;
  // Eigen::VectorXd beta = Eigen::VectorXd::Random(5);
  // Eigen::VectorXd intermed = X*beta;
  // res = stan::math::dot_product(beta,intermed);
  // std::vector<stan::math::var> beta_;

  // std::vector<double> g;
  // res.grad(beta,g);
  // Eigen::VectorXd analyticRes = (2.*(X*beta).array()).matrix();
  // std::cout << "autodiff = " << g.transpose() << std::endl;
  // std::cout << "analytic = " << analyticRes.transpose() << std::endl;

  // std::cout << stan::math::inv_gamma_log(2.,2.,2.) << std::endl;
  // int n = 5;
  // Eigen::MatrixXd X = Eigen::MatrixXd::Random(n,n);
  // X = stan::math::tcrossprod(X);

  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> X(5,5);
  X = Eigen::MatrixXd::Random(5,5);
  X = X.transpose() * X;
  quadratic f(X);
  Eigen::Matrix<double,Eigen::Dynamic,1> beta(5);
  for (int i = 0; i < 5; ++i) {
    beta(i) = 2.;
  }
  double fx;
  Eigen::Matrix<double,Eigen::Dynamic,1> grad_fx;
  stan::math::gradient(f,beta,fx,grad_fx);

  Eigen::Matrix<double, Eigen::Dynamic,1> analyticRes(5);
  // std::cout << f(beta) << std::endl;
  analyticRes = (2.*(X*beta).array()).matrix();
  std::cout << "autodiff = " << grad_fx.transpose() << std::endl;
  std::cout << "analytic = " << analyticRes.transpose() << std::endl;
  std::cout << "quad = " << f(beta) << std::endl;
  // int n = 5;
  // Eigen::VectorXd y(n);
  // Eigen::VectorXd mu(n);
  // Eigen::MatrixXd Sig(n,n);
  // for (int i = 0; i < n; ++i) {
  //   y(i) = 0.;
  //   mu(i) = 0.;
  //   for (int j = 0; j < n; ++j) {
  //     if (i==j) {
  //       Sig(i,i)= 3.;
  //     } else {
  //       Sig(i,j) = 0.;
  //     }
  //   }
  // }
  // std::cout << stan::math::multi_normal_log(y,mu,Sig) << std::endl;
  // int n = 5;
  // int J = 2;
  // int B = 5;
  // Eigen::MatrixXd diagm(B,B);
  // for (int i=0; i<B;++i) {
  //   for (int j=0;j<B;++j) {
  //     if (i==j) {
  //       diagm(i,i) = 1.;
  //     } else {
  //       diagm(i,j) = 0.;
  //     }
  //   }
  // }
  // // Eigen::MatrixXd parameters = ((Eigen::MatrixXd::Random(4+J,1).array()+1.)/2.).matrix();
  // Eigen::VectorXd parameters(4+J);
  // for (int k = 0; k < (4+J); ++k) {
  //   parameters(k) = 10.;
  // }
  // Eigen::VectorXd hyperparameters(7);
  // for (int k = 0; k < 7; ++k) {
  //   hyperparameters(k) = 10.;
  // }

  // // Eigen::MatrixXd hyperparameters = ((Eigen::MatrixXd::Random(7,1).array()+1.)/2.).matrix();
  // Eigen::MatrixXd y               = Eigen::MatrixXd::Zero(n,1);
  // // Eigen::MatrixXd mu              = Eigen::MatrixXd::Zero(n,1);
  // // Eigen::MatrixXd Sig             = Eigen::MatrixXd::Random(n,n);
  // // Sig = Sig.transpose() * Sig;
  // // std::cout << "log-multivariate normal="
  //           // << stan::math::multi_normal_log(y,mu,Sig)
  //           // << std::endl;
  // Eigen::MatrixXd vphi(n,J);
  // for (int i = 0; i < n; ++i) {
  //   for (int j = 0; j < J; ++j) {
  //     vphi(i,j) = 2.;
  //   }
  // }
  // // Eigen::MatrixXd vphi            = ((Eigen::MatrixXd::Random(n,J).array()+1.)/2.).matrix();
  // double tmp = logh<double>(parameters,hyperparameters,y,vphi,n,J,n);
  // // std::cout << "parameters=" << parameters.transpose() << std::endl;
  // // std::cout << "hyperparameters" << hyperparameters.transpose() << std::endl;
  // // std::cout << vphi << std::endl;
  // std::cout << "logh="
  //           << tmp
  //           << std::endl;
  // // std::cout << diagm << std::endl;
}