#include <iostream>
#include <cmath>
#include <armadillo>
#include "logH.h"

double logH(double p, double q, double r) {
  double mu0        = std::sqrt( ( p * r - 2.0 * q + std::sqrt( std::pow(( p * r - 2.0 * q ), 2.0) + 8.0 * q * r * ( p + 2.0 ) ) ) / ( 4.0 * q * r) );
  double sig0       = std::pow(-(-p / std::pow(mu0, 2.0) - 2.0 * q - 2.0 * ( 3.0 * r * std::pow(mu0, 2.0) + 1.0) / std::pow(( r * std::pow(mu0, 3.0) + mu0 ), 2.0)), -0.5);
  double hmu0       = p * std::log(mu0) - q * std::pow(mu0, 2.0) - std::log(r + std::pow(mu0, -2.0));
  double sig02      = sig0 + std::sqrt(2);
  double lowerLimit = -mu0 / sig02;
  double b          = 1;
  double epsilon    = 1;
  while(epsilon > 1.0e-5) {
    b *= 2;
    double maximum = std::max( std::exp( p * std::log(mu0 + b * sig02) - q * std::pow(mu0 + b * sig02, 2.0) - std::log(r + std::pow(mu0 + b * sig02, -2.0)) - hmu0 ), std::exp( p * std::log(mu0 - b * sig02) - q * std::pow(mu0 - b * sig02, 2.0) - std::log(r + std::pow(mu0 - b * sig02, -2.0)) - hmu0 ));
    double integrand_b = std::exp( p * std::log(mu0 + b * sig02) - q * std::pow(mu0 + b * sig02, 2.0) - std::log(r + std::pow(mu0 + b * sig02, -2.0)) - hmu0);
    epsilon = (-b > lowerLimit) ? maximum : integrand_b;
  }

  double s = 0.0;
  if (-b > lowerLimit) {
    double lower = -b;
    double upper = b;
    char   m     = 0;
    double h     = (upper - lower) / 200.0;
    double x_i, r_i;
    for (x_i = lower; x_i <= upper; x_i += h) {
      r_i = std::exp(p * std::log(mu0 + x_i * sig02) - q * std::pow(mu0 + x_i * sig02, 2.0) - std::log(r + std::pow(mu0 + x_i * sig02, -2.0)) - hmu0);
      if (x_i == lower || x_i == upper) {
        s += r_i;
      } else {
        m = !m;
        s += r_i * (m + 1) * 2.0;
      }
    }
    s *= h / 3.0;
  } else {
    double lower = lowerLimit;
    double upper = b;
    char   m     = 0;
    double h     = (upper - lower) / 100.0;
    double x_i, r_i;
    for (x_i = lower; x_i <= upper; x_i += h) {
      r_i = std::exp(p * std::log(mu0 + x_i * sig02) - q * std::pow(mu0 + x_i * sig02, 2.0) - std::log(r + std::pow(mu0 + x_i * sig02, -2.0)) - hmu0);
      if (x_i == lower || x_i == upper) {
        s += r_i;
      } else {
        m = !m;
        s += r_i * (m + 1) * 2.0;
      }
    }
    s *= h / 3.0;
  }
  return hmu0 + std::log(sig02) + std::log(s);
}
