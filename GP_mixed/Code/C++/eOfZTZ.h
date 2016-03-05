#ifndef _EOFZTZ_H
#define _EOFZTZ_H

#include <iostream>
#include <cmath>
#include <armadillo>

arma::mat eOfZTZ(arma::mat designMatrixX, arma::mat SMatrix, arma::colvec muOptimalLambda, arma::mat sigmaOptimalLambda);

#endif