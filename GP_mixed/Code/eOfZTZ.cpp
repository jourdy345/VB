#include <iostream>
#include <cmath>
#include <armadillo>
#include "eOfZTZ.h"

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