#include <iostream>
#include <cmath>
#include <armadillo>
#include "eOfZ.h"
#include "eOfZTZ.h"
#include "logH.h"
#include "lowerBound.h"

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