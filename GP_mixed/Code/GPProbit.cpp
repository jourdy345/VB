#include <iostream>
#include <armadillo>

arma::mat mvnRandomGenerate(int n, arma::vec mu, arma::mat sigma) {
   int nCols   = sigma.n_cols;
   arma::mat Y = arma::randn(n, nCols);

   return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

/* designMatrixX = X
   SMatrix       = S
   cutOff        = m
   dimension     = d
*/
arma::mat constructTmatrix(arma::mat designMatrixX, arma::mat SMatrix) {
  int nOfData         = designMatrixX.n_rows;
  int cutOff          = SMatrix.n_rows;
  int dimension       = SMatrix.n_cols;
  arma::mat Tmatrix   = arma::zeros<arma::mat>(nOfData * cutOff, dimension);

  for (int i = 1; i < nOfData + 1; i++) {
    Tmatrix.rows(cutOff * (i - 1), i * cutOff - 1) = SMatrix % arma::repmat(X.row(i - 1), cutOff, 1);
  }

  return Tmatrix;
}

/* nOfData = n
   cutOff  = m
*/
arma::mat eOfZ(arma::mat designMatrixX, arma::mat SMatrix, arma::vec muOptimalLambda, arma::mat sigmaOptimalLambda) {
  int nOfData = designMatrixX.n_rows;
  int cutOff  = SMatrix.n_rows;
  arma::mat Z = arma::zeros<arma::mat>(2 * cutOff, nOfData); 
  for (int i = 0; i < nOfData; i++) {
    arma::colvec ithColumn = Z.col(i);
    for (int r = 0; r < cutOff; r++) {
      arma::rowvec t_ir = SMatrix.row(r) % designMatrixX.row(i);
      ithColumn(r) = std::exp( -0.5 * arma::as_scalar( t_ir * sigmaOptimalLambda * t_ir.t() ) ) * std::cos( -0.5 * arma::as_scalar( t_ir * muOptimalLambda ) );
      ithColumn(r + cutOff) = std::exp( -0.5 * arma::as_scalar( t_ir * sigmaOptimalLambda * t_ir.t() ) ) * std::sin( -0.5 * arma::as_scalar( t_ir * muOptimalLambda ) );

    }
    Z.col(i) = ithColumn;
  }
  return Z;
}

arma::mat eOfZTZ(arma::mat designMatrixX, arma::mat SMatrix, arma::colvec muOptimalLambda, arma::mat sigmaOptimalLambda) {
  int nOfData      = designMatrixX.n_rows;
  int cutOff       = SMatrix.n_rows;
  arma::mat P      = arma::zeros<arma::mat>(cutOff, cutOff);
  arma::mat Q      = arma::zeros<arma::mat>(cutOff, cutOff);
  arma::mat R      = arma::zeros<arma::mat>(cutOff, cutOff);
  arma::mat P_temp = arma::zeros<arma::mat>(cutOff, cutOff);
  arma::mat Q_temp = arma::zeros<arma::mat>(cutOff, cutOff);
  arma::mat R_temp = arma::zeros<arma::mat>(cutOff, cutOff);

  for (int i = 0; i < nOfData; i++) {
    for (int r = 0; r < cutOff; r++) {
      for (int l = 0; l < cutOff; l++) {
        arma::rowvec tMinus_irl = (SMatrix.row(r) % designMatrixX.row(i)) - (SMatrix.row(l) % designMatrixX.row(i));
        arma::rowvec tPlus_irl  = (SMatrix.row(r) % designMatrixX.row(i)) + (SMatrix.row(l) % designMatrixX.row(i));
        P(r, l) = 0.5 * ( std::exp( -0.5 * arma::as_scalar( tMinus_irl * sigmaOptimalLambda * tMinus_irl.t() ) ) * std::cos( arma::as_scalar( tMinus_irl * muOptimalLambda ) ) + std::exp( arma::as_scalar( -0.5 * tPlus_irl * sigmaOptimalLambda * tPlus_irl.t() )) * std::cos( tPlus_irl * muOptimalLambda) );
        Q(r, l) = 0.5 * ( std::exp( -0.5 * arma::as_scalar( tMinus_irl * sigmaOptimalLambda * tMinus_irl.t() ) ) * std::sin( arma::as_scalar( tMinus_irl * muOptimalLambda ) ) + std::exp( arma::as_scalar( -0.5 * tPlus_irl * sigmaOptimalLambda * tPlus_irl.t() )) * std::sin( tPlus_irl * muOptimalLambda) );
        R(r, l) = 0.5 * ( std::exp( -0.5 * arma::as_scalar( tMinus_irl * sigmaOptimalLambda * tMinus_irl.t() ) ) * std::cos( arma::as_scalar( tMinus_irl * muOptimalLambda ) ) - std::exp( arma::as_scalar( -0.5 * tPlus_irl * sigmaOptimalLambda * tPlus_irl.t() )) * std::cos( tPlus_irl * muOptimalLambda) );
      }
    }
    P_temp += P;
    Q_temp += Q;
    R_temp += R;
  }
  arma::mat finalMatrix   = arma::join_rows(P, Q);
  arma::mat finalMatrix_2 = arma::join_rows(Q.t(), R);
  arma::mat final         = arma::join_cols(finalMatrix, finalMatrix_2);

  return final;
}

// Numerical integration of logH (the integration function was replaced with Simpson's rule)
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


double lowerBound(arma::colvec yResponse, arma::mat designMatrixX, arma::mat SMatrix, arma::colvec muOptimalLambda, arma::mat sigmaOptimalLambda, arma::colvec muOptimalAlpha, arma::mat sigmaOptimalAlpha, arma::colvec muOptimalBeta, arma::mat sigmaOptimalBeta, arma::mat A) {
  double part1 = -0.5 * (arma::trace( eOfZTZ( designMatrixX, SMatrix, muOptimalLambda, sigmaOptimalLambda ) * sigmaOptimalAlpha ) + arma::as_scalar(muOptimalAlpha.t() * ( eOfZTZ( designMatrixX, SMatrix, muOptimalLambda, sigmaOptimalLambda) - (eOfZ(designMatrixX, SMatrix, muOptimalLambda, sigmaOptimalLambda).t() * eOfZ(designMatrixX, SMatrix, muOptimalLambda, sigmaOptimalLambda)) ) * muOptimalAlpha) + arma::trace(A.t() * A  * sigmaOptimalBeta));
}

int main() {
  std::cout << logH(20, 10, 10);

  return 0;
}
