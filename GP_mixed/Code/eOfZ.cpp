#include <iostream>
#include <cmath>
#include <armadillo>
#include "eOfZ.h"

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
  return Z.t();
}
