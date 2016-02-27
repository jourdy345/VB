#ifndef _LOWERBOUND_H
#define _LOWERBOUND_H

#include <iostream>
#include <cmath>
#include <armadillo>

double lowerBound(arma::colvec yResponse, arma::mat designMatrixX, arma::mat SMatrix, arma::colvec muOptimalLambda, arma::mat sigmaOptimalLambda, arma::colvec muOptimalAlpha, arma::mat sigmaOptimalAlpha, arma::colvec muOptimalBeta, arma::mat sigmaOptimalBeta, arma::mat A, arma::colvec priorMuBeta, arma::mat priorSigmaBeta, arma::colvec priorMuLambda, arma::mat priorSigmaLambda, double A_sigma, double C_sigma);

#endif