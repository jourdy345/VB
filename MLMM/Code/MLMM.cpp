// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <Rmath.h>
using namespace Rcpp; using namespace arma;

// [[Rcpp::export]]
double median_rcpp(NumericVector x) {
   NumericVector y = clone(x);
   int n, half;
   double y1, y2;
   n = y.size();
   half = n / 2;
   if(n % 2 == 1) {
      // median for odd length vector
      std::nth_element(y.begin(), y.begin()+half, y.end());
      return y[half];
   } else {
      // median for even length vector
      std::nth_element(y.begin(), y.begin()+half, y.end());
      y1 = y[half];
      std::nth_element(y.begin(), y.begin()+half-1, y.begin()+half);
      y2 = y[half-1];
      return (y1 + y2) / 2.0;
   }
}

// [[Rcpp::export]]
double fit_t (NumericVector x) {
  int n = x.size();
  double v_old = median_rcpp(x);
  double v_new;
  int step = 0;
  int step2 = 0;
  while (TRUE) {
    step += 1;
    double l_first_der = 0;
    double l_second_der = 0;
    for (int i = 0; i < n; ++i) {
      l_first_der += R::digamma((v_old+1) / 2) - 1 / (2 * v_old) - R::digamma(v_old / 2) - 0.5 * log(1 + (x(i) * x(i)) / v_old) + 0.5 * (v_old + 1) * (pow(x(i), 2) / (pow(v_old, 2) + v_old * pow(x(i), 2)));
      l_second_der += R::trigamma((v_old+1) / 2) + 1 / (2 * pow(v_old, 2)) - R::trigamma(v_old / 2) + (pow(x(i), 2) / (pow(v_old, 2) + v_old * pow(x(i), 2))) - 0.5 * (v_old + 1) * pow(x(i), 2) * (2 * v_old + pow(x(i), 2)) / pow((pow(v_old, 2) + v_old * pow(x(i), 2)), 2);
    }
    v_new = v_old - (l_first_der / l_second_der);
    if (std::abs(v_new - v_old) < 0.0001) {
      break;
    } else if (step > 50) {
      step = 0;
      v_old = R::mean(x);
    } else {
      v_old = v_new;
    }
  }
  return v_new;
}

// [[Rcpp::export]]
double lower_bound (Rcpp::List Sigma_beta, Rcpp:List Sigma_beta_q, arma::mat mu_beta_q, Rcpp::List Sigma_b_q, arma::mat mu_b_q, arma::vec alpha_b, arma::vec lambda_b, arma::vec alpha_b_q, arma::vec lambda_b_q, arma::vec alpha_a, arma::vec lambda_a, arma::vec alpha_a_q, arma::vec lambda_a_q, arma::mat q_ij, ) {

}



