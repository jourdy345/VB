#ifndef _MVNRANDOMGENERATE_H
#define _MVNRANDOMGENERATE_H

#include <iostream>
#include <cmath>
#include <armadillo>

arma::mat mvnRandomGenerate(int n, arma::vec mu, arma::mat sigma);

#endif