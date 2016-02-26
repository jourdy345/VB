// don't forget to add armadillo library to linker
// g++ -c GPProbit.cpp
// g++ -o GPProbit GPProbit.o -llapack -lblas -larmadillo

#define _USE_MATH_DEFINES
#include <iostream>
#include <armadillo>
#include <cmath> //for erfc function and math constants

arma::mat mvnRandomGenerate(int n, arma::vec mu, arma::mat sigma) {
   int nCols   = sigma.n_cols;
   arma::mat Y = arma::randn(n, nCols);

   return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
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
      arma::colvec t_ir = (SMatrix.row(r) % designMatrixX.row(i)).t();
      ithColumn(r) = std::exp( -0.5 * arma::dot( t_ir, sigmaOptimalLambda * t_ir ) ) * cos( -0.5 * arma::dot( t_ir, muOptimalLambda ) );
      ithColumn(r + cutOff) = std::exp( -0.5 * arma::dot( t_ir, sigmaOptimalLambda * t_ir ) ) * sin( -0.5 * arma::dot( t_ir, muOptimalLambda ) );

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
        arma::colvec tMinus_irl = ((SMatrix.row(r) % designMatrixX.row(i)) - (SMatrix.row(l) % designMatrixX.row(i))).t();
        arma::colvec tPlus_irl  = ((SMatrix.row(r) % designMatrixX.row(i)) + (SMatrix.row(l) % designMatrixX.row(i))).t();
        P(r, l) = 0.5 * ( std::exp( -0.5 * ( arma::dot(tMinus_irl, sigmaOptimalLambda * tMinus_irl) )) * cos( arma::dot(tMinus_irl, muOptimalLambda) ) + std::exp( -0.5 * arma::dot(tPlus_irl, sigmaOptimalLambda * tPlus_irl )) * cos( arma::dot(tPlus_irl, muOptimalLambda) ));
        Q(r, l) = 0.5 * ( std::exp( -0.5 * ( arma::dot(tMinus_irl, sigmaOptimalLambda * tMinus_irl) )) * sin( arma::dot(tMinus_irl, muOptimalLambda) ) + std::exp( -0.5 * arma::dot(tPlus_irl, sigmaOptimalLambda * tPlus_irl )) * sin( arma::dot(tPlus_irl, muOptimalLambda) ));
        R(r, l) = 0.5 * ( std::exp( -0.5 * ( arma::dot(tMinus_irl, sigmaOptimalLambda * tMinus_irl) )) * cos( arma::dot(tMinus_irl, muOptimalLambda) ) - std::exp( -0.5 * arma::dot(tPlus_irl, sigmaOptimalLambda * tPlus_irl )) * cos( arma::dot(tPlus_irl, muOptimalLambda) ));
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


double lowerBound(arma::colvec yResponse, arma::mat designMatrixX, arma::mat SMatrix, arma::colvec muOptimalLambda, arma::mat sigmaOptimalLambda, arma::colvec muOptimalAlpha, arma::mat sigmaOptimalAlpha, arma::colvec muOptimalBeta, arma::mat sigmaOptimalBeta, arma::mat A, arma::colvec priorMuBeta, arma::mat priorSigmaBeta, arma::colvec priorMuLambda, arma::mat priorSigmaLambda, double A_sigma, double C_sigma) {
  int nOfData   = yResponse.size();
  int cutOff    = SMatrix.n_rows;
  int dimension = SMatrix.n_cols;

  double n      = static_cast<double>(nOfData);
  double m      = static_cast<double>(cutOff);
  double d      = static_cast<double>(dimension);
  double s      = A.n_cols;

  double part1 = -0.5 * (arma::trace( eOfZTZ( designMatrixX, SMatrix, muOptimalLambda, sigmaOptimalLambda ) * sigmaOptimalAlpha ) + arma::as_scalar(muOptimalAlpha.t() * ( eOfZTZ( designMatrixX, SMatrix, muOptimalLambda, sigmaOptimalLambda) - (eOfZ(designMatrixX, SMatrix, muOptimalLambda, sigmaOptimalLambda).t() * eOfZ(designMatrixX, SMatrix, muOptimalLambda, sigmaOptimalLambda)) ) * muOptimalAlpha) + arma::trace(A.t() * A  * sigmaOptimalBeta));

  double part2 = 0.0;
  arma::colvec functand = eOfZ(designMatrixX, SMatrix, muOptimalLambda, sigmaOptimalLambda) * muOptimalAlpha + A * muOptimalBeta;
  for (int index = 0; index < nOfData; index++) {
    if (yResponse(index) == 1.0) {
      part2 += std::log(0.5 * std::erfc(-functand(index) * M_SQRT1_2));
    } else {
      part2 += std::log( 1 - 0.5 * std::erfc(-functand(index) * M_SQRT1_2) );
    }
  }
  double sign, detSigmaOptimalAlpha, detSigmaOptimalBeta, detSigmaOptimalLambda, detPriorSigmaLambda, detPriorSigmaBeta;
  arma::log_det(detSigmaOptimalAlpha, sign, sigmaOptimalAlpha);
  arma::log_det(detSigmaOptimalBeta, sign, sigmaOptimalBeta);
  arma::log_det(detSigmaOptimalLambda, sign, sigmaOptimalLambda);
  arma::log_det(detPriorSigmaLambda, sign, priorSigmaLambda);
  arma::log_det(detPriorSigmaBeta, sign, priorSigmaBeta);

  double part3 = m * std::log(m) + 0.5 * detSigmaOptimalAlpha + m - 0.5 * detPriorSigmaBeta - 0.5 * ( arma::trace(arma::solve(priorSigmaBeta, sigmaOptimalBeta)) + arma::dot((muOptimalBeta - priorMuBeta), arma::solve(priorSigmaBeta, muOptimalBeta - priorMuBeta)));
  double part4 = 0.5 * (detSigmaOptimalBeta + s) + std::log(2.0 * A_sigma) - std::log(M_PI) + logH(2.0 * m - 2.0, C_sigma, std::pow(A_sigma, 2.0)) + 0.5 * (detSigmaOptimalLambda + d);
  double part5 = -0.5 * ( detPriorSigmaLambda + arma::trace(arma::solve(priorSigmaLambda, sigmaOptimalLambda)) + arma::dot((muOptimalLambda - priorMuLambda), arma::solve(priorSigmaLambda, muOptimalLambda - priorMuLambda)));

  return part1 + part2 + part3 + part4 + part5;
}

void sparseGPVBProbit(arma::colvec yResponse, arma::mat designMatrixX, arma::mat SMatrix, arma::mat A, arma::colvec priorMuBeta, arma::mat priorSigmaBeta, arma::colvec priorMuLambda, arma::mat priorSigmaLambda, double A_sigma) {
  int nOfData   = yResponse.size();
  int cutOff    = SMatrix.n_rows;
  int dimension = SMatrix.n_cols;
  int sInt      = A.n_cols;
  double n      = static_cast<double>(nOfData);
  double m      = static_cast<double>(cutOff);
  double d      = static_cast<double>(dimension);
  double s      = static_cast<double>(sInt);
  double oldLowerBound = -10000.0;
  double newLowerBound;
  // initialize variational parameters
  arma::mat sigmaOptimalLambda; sigmaOptimalLambda.eye(dimension, dimension);
  arma::colvec muOptimalLambda; muOptimalLambda.randn(dimension);
  arma::colvec muOptimalBeta; muOptimalBeta.randn(sInt);
  arma::mat sigmaOptimalAlpha, sigmaOptimalBeta;
  arma::colvec muOptimalResponseStar, muOptimalAlpha;
  double C_sigma = 10;
  while (true) {
    // update μ(α), Σ(α)
    arma::mat tempDiagonalMatrix(2 * cutOff, 2 * cutOff); tempDiagonalMatrix.eye();
    sigmaOptimalAlpha = ( m * std::exp(logH(2.0 * m - 4.0, C_sigma, std::pow(A_sigma, 2.0)) - logH(2.0 * m - 2.0, C_sigma, std::pow(A_sigma, 2.0))) * tempDiagonalMatrix + eOfZTZ( designMatrixX, SMatrix, muOptimalLambda, sigmaOptimalLambda ) ).i();
    muOptimalAlpha    = sigmaOptimalAlpha * ( eOfZ( designMatrixX, SMatrix, muOptimalLambda, sigmaOptimalLambda ).t() * muOptimalResponseStar - ( A * eOfZ( designMatrixX, SMatrix, muOptimalLambda, sigmaOptimalLambda )).t() * muOptimalBeta );

    // update μ(y*)
    arma::colvec functand = eOfZ( designMatrixX, SMatrix, muOptimalLambda, sigmaOptimalLambda ) * muOptimalAlpha + A * muOptimalBeta;
    for (int index = 0; index < nOfData; index++) {
      if (yResponse(index) == 1.0) {
        muOptimalResponseStar(index) = functand(index) + ((1.0 / std::sqrt(2.0 * M_PI)) * std::exp(-0.5 * functand(index) * functand(index))) / (0.5 * std::erfc(-functand(index) * M_SQRT1_2));
      } else {
        muOptimalResponseStar(index) = functand(index) + ((1.0 / std::sqrt(2.0 * M_PI)) * std::exp(-0.5 * functand(index) * functand(index))) / (0.5 * std::erfc(-functand(index) * M_SQRT1_2) - 1.0);
      }
    }

    // update μ(β), Σ(β)
    sigmaOptimalBeta = (priorSigmaBeta.i() + A.t() * A).i();
    muOptimalBeta = sigmaOptimalBeta * ( arma::solve( priorSigmaBeta, priorMuBeta ) + A.t() * muOptimalResponseStar - A.t() * eOfZ( designMatrixX, SMatrix, muOptimalLambda, sigmaOptimalLambda ) * muOptimalAlpha );

    // update C(σ)
    C_sigma = 0.5 * m * (arma::trace(sigmaOptimalAlpha) + arma::dot(muOptimalAlpha, muOptimalAlpha));

    // update μ(λ), Σ(λ)
    arma::mat F1(dimension, dimension);
    arma::mat F2(dimension, dimension);
    arma::colvec F3(dimension);
    arma::colvec F4(dimension);
    arma::colvec t_ij(dimension);
    arma::colvec tPlus_ijr(dimension);
    arma::colvec tMinus_ijr(dimension);
    arma::colvec AmuOptimalBeta = A * muOptimalBeta;
    arma::mat ABD = sigmaOptimalAlpha + muOptimalAlpha * muOptimalAlpha.t();
    arma::mat A = ABD.submat(0, 0, cutOff - 1, cutOff - 1);
    arma::mat B = ABD.submat(cutOff - 1, 0, 2 * cutOff - 1, cutOff - 1);
    arma::mat D = ABD.submat(cutOff - 1, cutOff - 1, 2 * cutOff - 1, 2 * cutOff - 1);
    for (int i = 0; i < nOfData; i++) {
      double temp = muOptimalResponseStar(i) - AmuOptimalBeta(i);
      for (int j = 0; j < cutOff; j++) {
        t_ij = (SMatrix.row(j) % designMatrixX.row(i)).t();
        F1 = temp * std::exp( -0.5 * arma::dot(t_ij, sigmaOptimalLambda * t_ij) ) * (muOptimalAlpha(j) * cos( arma::dot(t_ij, muOptimalLambda) ) + muOptimalAlpha(j + cutOff) * sin( arma::dot( t_ij, muOptimalLambda ) )) * (t_ij * t_ij.t());
        F3 = temp * std::exp( -0.5 * arma::dot(t_ij, sigmaOptimalLambda * t_ij) ) * (muOptimalAlpha(j + cutOff) * cos( arma::dot(t_ij, muOptimalLambda) ) - muOptimalLambda(j) * sin( arma::dot( t_ij, muOptimalLambda ) )) * t_ij;
        for (int r = 0; r < cutOff; r++) {
          tPlus_ijr = (SMatrix.row(j) % designMatrixX.row(i) + SMatrix.row(r) % designMatrixX.row(i)).t();
          tMinus_ijr = (SMatrix.row(j) % designMatrixX.row(i) - SMatrix.row(r) % designMatrixX.row(i)).t();
          F2 = (std::exp( -0.5 * arma::dot( tMinus_ijr, sigmaOptimalLambda * tMinus_ijr ) ) * ( (A(j, r) + D(j, r) ) * cos( arma::dot( tMinus_ijr, muOptimalLambda ) ) + 2.0 * B(j, r) * sin( arma::dot( tMinus_ijr, muOptimalLambda ) ) ) * (tMinus_ijr * tMinus_ijr.t())) + (std::exp( -0.5 * arma::dot( tPlus_ijr, sigmaOptimalLambda * tPlus_ijr ) ) * ( (A(j, r) - D(j, r)) * cos( arma::dot(tPlus_ijr, muOptimalLambda) ) + 2.0 * B(j, r) * sin( arma::dot(tPlus_ijr, muOptimalLambda) ) ) * (tPlus_ijr * tPlus_ijr.t()));
          F4 = (std::exp( -0.5 * arma::dot( tMinus_ijr, sigmaOptimalLambda * tMinus_ijr ) ) * ( 2.0 * B(j, r) * cos( arma::dot(tMinus_ijr, muOptimalLambda) ) - (A(j, r) + D(j, r)) * sin( arma::dot( tMinus_ijr, muOptimalLambda ) ) ) * tMinus_ijr) + std::exp( -0.5 * arma::dot( tPlus_ijr, sigmaOptimalLambda * tPlus_ijr ) ) * ( 2.0 * B(j, r) * cos( arma::dot(tPlus_ijr, muOptimalLambda) ) + (D(j, r) - A(j, r)) * sin( arma::dot( tPlus_ijr, muOptimalLambda ) ) ) * tPlus_ijr;
        }
      }
    }
    F2 *= -0.5;
    F4 *= -0.25;

    sigmaOptimalLambda = (priorSigmaLambda.i() + F1 + F2).i();
    muOptimalLambda   += sigmaOptimalLambda * (arma::solve(priorSigmaLambda, priorMuLambda - muOptimalLambda) + F3 + F4);

    newLowerBound = lowerBound(yResponse, designMatrixX, SMatrix, muOptimalLambda, sigmaOptimalLambda, muOptimalAlpha, sigmaOptimalAlpha, muOptimalBeta, sigmaOptimalBeta, A, priorMuBeta, priorSigmaBeta, priorMuLambda, priorSigmaLambda, A_sigma, C_sigma);
    if (newLowerBound - oldLowerBound < 0e-7) break;
    else oldLowerBound = newLowerBound;
  }
  std::cout << "μ(α)" << std::endl;
  std::cout << muOptimalAlpha.t() << std::endl;
  std::cout << "Σ(α)" << std::endl;
  std::cout << sigmaOptimalAlpha << std::endl;
  std::cout << "μ(β)" << std::endl;
  std::cout << muOptimalBeta.t() << std::endl;
  std::cout << "Σ(β)" << std::endl;
  std::cout << sigmaOptimalBeta << std::endl;
  std::cout << "μ(λ)" << std::endl;
  std::cout << muOptimalLambda.t() << std::endl;
  std::cout << "Σ(λ)" << std::endl;
  std::cout << sigmaOptimalLambda << std::endl;
}

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
