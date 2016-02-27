#ifndef _SPARSEGPVBPROBIT_H
#define _SPARSEGPVBPROBIT_H


void sparseGPVBProbit(arma::colvec yResponse, arma::mat designMatrixX, arma::mat SMatrix, arma::mat A, arma::colvec priorMuBeta, arma::mat priorSigmaBeta, arma::colvec priorMuLambda, arma::mat priorSigmaLambda, double A_sigma);

#endif