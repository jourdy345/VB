#ifndef _EOFZ_H
#define _EOFZ_H

#include <iostream>
#include <cmath>
#include <armadillo>

arma::mat eOfZ(arma::mat designMatrixX, arma::mat SMatrix, arma::vec muOptimalLambda, arma::mat sigmaOptimalLambda);

#endif