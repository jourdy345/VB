#define _USE_MATH_DEFINES
#include "./math/stan/math.hpp"
#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>
using Eigen::Matrix;
using Eigen::Dynamic;

/* sudo g++ -I ~/Desktop/tmp/VB_cpp/math -I ~/Desktop/tmp/VB_cpp/math/lib/eigen_3.2.9/ -I ~/Desktop/tmp/VB_cpp/math/lib/boost_1.62.0/ -I ~/Desktop/tmp/VB_cpp/math/lib/cvodes_2.9.0/include/ -std=c++11 -stdlib=libc++ VB_reg_f1.cpp -o VB_reg_f1.out */


// function that outputs data to csv
const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
void writeToCSVfile(std::string name, Eigen::MatrixXd matrix)
{
    std::ofstream file(name.c_str());
    file << matrix.format(CSVFormat);
 }


// double logh_fun(Matrix<double,Dynamic,1>& parameters, Matrix<double,Dynamic,1>& hyperparameters,
//                 Matrix<double,Dynamic,1>& y, Matrix<double,Dynamic,Dynamic>& vphi,
//                 int& n, int& J)
// {
//   Matrix<double,Dynamic,1> Diagm(n,n);
//   for (int i = 0; i < n; ++i) {
//     for (int j = 0; j < n; ++j) {
//       if (i == j) {
//         Diagm(i,i) = 1.;
//       } else {
//         Diagm(i,j) = 0.;
//       }
//     }
//   }
//   Matrix<double,Dynamic,1> theta(J);
//   for (int j = 0; j < J; ++j) {
//     theta(j) = parameters(j);
//   }
//   double psi  = parameters(J);
//   double sig2 = parameters(J+1);
//   double tau2 = parameters(J+2);
//   double gam  = parameters(J+3);

//   // extract hyperparameters
//   double r_sig = hyperparameters(0);
//   double s_sig = hyperparameters(1);
//   double r_tau = hyperparameters(2);
//   double s_tau = hyperparameters(3);
//   double     a = hyperparameters(4);
//   double     b = hyperparameters(5);
//   double    w0 = hyperparameters(6);

//   double sigma2;
//   Matrix<double,Dynamic,1> mu(J+4);
//   Matrix<double,Dynamic,Dynamic> Sig(J+4,J+4);
//   // Calculate log-likelihood+log(prior)
//   mu = vphi * theta;
//   Sig = (sig2*Diagm.array()).matrix();
//   double res = stan::math::multi_normal_log(y,mu,Sig)+stan::math::inv_gamma_log(sig2,r_sig/2.,s_sig/2.)+stan::math::inv_gamma_log(tau2,r_tau/2.,s_tau/2.)+stan::math::inv_gamma_log(psi,a,b)+std::log(w0)-w0*gam;
//   for (int j = 0; j < J; ++j) {
//     sigma2 = sig2*tau2*stan::math::exp(-(j+1.)*gam);
//     res += stan::math::normal_log(theta(j),0.,sigma2);
//   }
//   return res;
// }


struct logh {
  const Matrix<double,Dynamic,1> hyperparameters_;
  const Matrix<double,Dynamic,1> y_;
  const Matrix<double,Dynamic,Dynamic> vphi_;
  const Matrix<double,Dynamic,Dynamic> Diagm_;
  const int n_;
  const int J_;
  logh(const Matrix<double,Dynamic,1>& hyperparameters,
       const Matrix<double,Dynamic,1>& y,
       const Matrix<double,Dynamic,Dynamic>& vphi,
       const Matrix<double,Dynamic,Dynamic>& Diagm,
       const int& n,
       const int& J) : hyperparameters_(hyperparameters),
                       y_(y),
                       vphi_(vphi),
                       Diagm_(Diagm),
                       n_(n),
                       J_(J) { }
  template <typename T>
  T operator()(const Matrix<T,Dynamic,1>& parameters) const {
    // std::cout << "START" << std::endl;
    // Matrix<T,Dynamic,1> Diagm(n_,n_);
    // Diagm.setIdentity();
    // std::cout << "START2" << std::endl;
    Matrix<T,Dynamic,1> theta(J_);
    for (int j = 0; j < J_; ++j) {
      theta(j) = parameters(j);
    }
    T psi  = parameters(J_);
    T sig2 = parameters(J_+1);
    T tau2 = parameters(J_+2);
    T gam  = parameters(J_+3);

    // extract hyperparameters
    T r_sig = hyperparameters_(0);
    T s_sig = hyperparameters_(1);
    T r_tau = hyperparameters_(2);
    T s_tau = hyperparameters_(3);
    T     a = hyperparameters_(4);
    T     b = hyperparameters_(5);
    T    w0 = hyperparameters_(6);

    T sigma2;
    Matrix<T,Dynamic,1> mu(J_+4);
    Matrix<T,Dynamic,Dynamic> Sig(J_+4,J_+4);
    // Calculate log-likelihood+log(prior)
    mu = stan::math::multiply(vphi_,theta);
    // Sig = stan::math::multiply(sig2,Diagm_);
    T   res = stan::math::inv_gamma_log(sig2,r_sig/2.,s_sig/2.)+stan::math::inv_gamma_log(tau2,r_tau/2.,s_tau/2.)+stan::math::inv_gamma_log(psi,a,b)+stan::math::log(w0)-w0*gam;
    for (int i = 0; i < n_; ++i) {
      res += stan::math::normal_log(y_(i),mu(i),sig2);
    }
    for (int j = 0; j < J_; ++j) {
      sigma2 = sig2*tau2*stan::math::exp(-(j+1.)*gam);
      res += stan::math::normal_log(theta(j),0.,sigma2);
    }
    return res;
  }
};




int main() {
  int n = 200;
  int J = 35;
  // set RNG
  std::random_device rd;
  std::mt19937 gen(rd());

  std::normal_distribution<double> randn(0.,1.);
  std::uniform_real_distribution<double> rand(0.,1.);

  Matrix<double,Dynamic,1> x(n);
  for (int i = 0; i < n; ++i) {
    x(i) = rand(gen);
  }

  Matrix<double,Dynamic,1> fx(n);
  for (int i = 0; i < n; ++i) {
    fx(i) = 5.+stan::math::exp(10.*x(i)-5.)/(1.+std::exp(10.*x(i)-5.));
  }

  double sig2_true = 3.4;
  double xmin = stan::math::min(x);
  double xmax = stan::math::max(x);
  x = stan::math::elt_divide(stan::math::subtract(x,xmin),xmax-xmin);

  // response variable
  Matrix<double,Dynamic,1> y(n);
  for (int i = 0; i < n; ++i) {
    y(i) = fx(i)+randn(gen)*stan::math::sqrt(sig2_true);
  }

  Matrix<double,1,Dynamic> linspace(J);
  for (int j = 0; j < J; ++j) {
    linspace(j) = M_PI*double(j+1);
  }
  Matrix<double,Dynamic,Dynamic> vphi = ((x*linspace).array().cos() * std::sqrt(2.)).matrix();

  Matrix<double,Dynamic,1> hyperparameters(7);
  for (int i = 0; i < 6; ++i) {
    hyperparameters(i) = 10.;
  }
  hyperparameters(6) = 2.;


  // set variational parameters & auxiliary variables
  Matrix<double,Dynamic,1> s_k(J+4);
  Matrix<double,Dynamic,1> rho(J+4);
  Matrix<double,Dynamic,1> mu_post(J+4);
  Matrix<double,Dynamic,Dynamic> L_post(J+4,J+4);
  for (int i = 0; i < (J+4); ++i) {
    mu_post(i) = 0.;
        s_k(i) = 1.;
        rho(i) = 0.;
    for (int j = 0; j < (J+4); ++j) {
      if (i == j) {
        L_post(i,i) = 1.;
      } else {
        L_post(i,j) = 0.;
      }
    }
  }

  double     psi;
  double     sig2;
  double     tau2;
  double     gam;
  double     logh_v;


  double      LB = 0.;
  int          S = 20;
  double     d_S = double(S);
  double epsilon = 1.0e-16;
  double     alp = 0.3;
  double       t = 1.;
  bool      cont = true;
  int    iterMax = 10000;
  int      nIter = 0;
  int      count = 0;
  int    tmp_int = 0;
  double     tmp = 0.;
  double tmp_pow = -0.5+epsilon;
  double     diff = 0.;
  std::vector<double> LBC;
  std::vector<double>::iterator Iter;

  Matrix<double,Dynamic,Dynamic> thetaContainer(J+4,iterMax);
  Matrix<double,Dynamic,Dynamic> L_1t(J+4,J+4);
  Matrix<double,Dynamic,Dynamic> z_eta(S,J+4);
  Matrix<double,Dynamic,Dynamic> Sig(J+4,J+4);
  Matrix<double,Dynamic,Dynamic> Diag(n,n);
  Matrix<double,Dynamic,1>       theta(J);
  Matrix<double,Dynamic,1>       z(J+4);
  Matrix<double,Dynamic,1>       zeta(J+4);
  Matrix<double,Dynamic,1>       Theta(J+4);
  Matrix<double,Dynamic,1>       grad_Theta;
  Matrix<double,Dynamic,1>       nabla_Tinv(J+4);
  Matrix<double,Dynamic,1>       nabla_logdet(J+4);
  Matrix<double,Dynamic,1>       grad_tmp(J+4);
  nabla_Tinv.fill(1.);
  nabla_logdet.setZero();
  Diag.setIdentity();

  logh log_likelihood_prior(hyperparameters,y,vphi,Diag,n,J);
  Theta.fill(2.);
  std::cout << log_likelihood_prior(Theta) << std::endl;
  // initialize log_likelihood_prior as a functor
  // tmp = log_likelihood_prior(Theta);
  while (cont) {
    for (int i = 0; i < S; ++i) {
      for (int j = 0; j < (J+4); ++j) {
        z_eta(i,j) = randn(gen);
      }
    }

    Matrix<double,Dynamic,1>       grad_mu(J+4);
    Matrix<double,Dynamic,Dynamic> grad_L(J+4,J+4);
    grad_mu.setZero();
    grad_L.setZero();
    for (int s = 0; s < S; ++s) {
      z = z_eta.row(s).transpose();
      zeta = L_post*z+mu_post;
      for (int j = 0; j < J; ++j) {
        theta(j) = zeta(j);
      }
      psi  = std::log(std::exp(zeta(J))+1.);
      sig2 = std::log(std::exp(zeta(J+1))+1.);
      tau2 = std::log(std::exp(zeta(J+2))+1.);
      gam  = std::log(std::exp(zeta(J+3))+1.);
      for (int i = 0; i < J; ++i) {
        Theta(i) = theta(i);
      }
      Theta(J)   = psi;
      Theta(J+1) = sig2;
      Theta(J+2) = tau2;
      Theta(J+3) = gam;
      stan::math::gradient(log_likelihood_prior,Theta,logh_v,grad_Theta);
      for (int i = 0; i < 4; ++i) {
        nabla_Tinv(J+i)   = std::exp(zeta(J+i))/(1.+std::exp(zeta(J+i)));
        nabla_logdet(J+i) = 1./(1.+std::exp(zeta(J+i)));
      }
      grad_tmp = ((grad_Theta.array()/d_S) * nabla_Tinv.array()).matrix()+stan::math::elt_divide(nabla_logdet,d_S);
      grad_mu += grad_tmp;
      L_1t     = stan::math::inverse(L_post).transpose();
      grad_L  += grad_tmp * z.transpose()+(L_1t.array()/d_S).matrix();
    }
    if (nIter == 0) {
      s_k = (grad_mu.array()*grad_mu.array()).matrix();
    } else {
      s_k = ((alp*grad_mu.array())*grad_mu.array()).matrix()+((1.-alp)*s_k.array()).matrix();
    }
    tmp_int = nIter+1;
    tmp = double(tmp_int);
    rho = ((0.1*std::pow(tmp,tmp_pow)/(t+s_k.array().sqrt())).matrix());
    mu_post += (rho.array()*grad_mu.array()).matrix();
    L_post  += stan::math::diag_matrix(rho)*grad_L;
    LB = 0.;

    for (int s = 0; s < S; ++s) {
      tmp = 0.;
      for (int i = 0; i < (J+4); ++i) {
        z(i) = randn(gen);
      }
      zeta = L_post*z+mu_post;
      for (int j = 0; j < J; ++j) {
        theta(j) = zeta(j);
      }
      psi  = std::log(std::exp(zeta(J))+1.);
      sig2 = std::log(std::exp(zeta(J+1))+1.);
      tau2 = std::log(std::exp(zeta(J+2))+1.);
      gam  = std::log(std::exp(zeta(J+3))+1.);
      for (int i = 0; i < J; ++i) {
        Theta(i) = theta(i);
      }
      Theta(J)   = psi;
      Theta(J+1) = sig2;
      Theta(J+2) = tau2;
      Theta(J+3) = gam;

      Sig = stan::math::tcrossprod(L_post);
      for (int j = J; j < (J+4); ++j) {
        tmp += zeta(j)-stan::math::log_sum_exp(0.,zeta(j));
      }
      LB += log_likelihood_prior(Theta)/d_S;
      LB += (tmp-stan::math::multi_normal_log(zeta,mu_post,Sig))/d_S;
    }
    std::cout << nIter+1 << "th LB=" << LB << std::endl;
    LBC.push_back(LB);
    thetaContainer.col(nIter) = mu_post;
    if ( (nIter > 200) && ((nIter+1) % 100 == 0) ) {
      diff = 0.;
      for (Iter = LBC.end(); Iter != (LBC.end()-100); --Iter) {
        diff += (*Iter - *(Iter-100));
      }
      diff /= 100.;
      if ((diff < 0.) || (stan::math::abs(diff) < 0.001)) {
        count += 1;
        if (count > 3) {
          cont = false;
        }
      } else {
        count = 0;
      }
    }
    nIter += 1;
    if (nIter == iterMax) {
      cont = false;
    }
  }

  
  writeToCSVfile("y_f1",y);
  writeToCSVfile("fx_f1",fx);
  writeToCSVfile("x_f1",x);
  writeToCSVfile("vphi_f1",vphi);
  writeToCSVfile("thetaContainer_f1",thetaContainer);
  writeToCSVfile("L_post_f1",L_post);

  return 0;
}