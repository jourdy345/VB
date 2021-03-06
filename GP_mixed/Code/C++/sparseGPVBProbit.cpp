
#include <iostream>
#include <armadillo>
#include <cmath>
#include "eOfZ.h"
#include "eOfZTZ.h"
#include "logH.h"
#include "lowerBound.h"
#include "sparseGPVBProbit.h"

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
  double newLowerBound = 0.0;
  double dif = 100.0;
  // initialize variational parameters
  arma::mat sigmaOptimalLambda; sigmaOptimalLambda.eye(dimension, dimension);
  arma::colvec muOptimalLambda; muOptimalLambda.randn(dimension);
  arma::colvec muOptimalBeta; muOptimalBeta.randn(sInt);
  arma::colvec muOptimalAlpha; muOptimalAlpha.randn(2 * cutOff);
  arma::mat sigmaOptimalAlpha(2 * cutOff, 2 * cutOff);
  arma::mat sigmaOptimalBeta(sInt, sInt);
  arma::colvec muOptimalResponseStar(nOfData);
  double C_sigma = 10.0;
  int count = 0;
  int count2 = 0;
  while ((dif > 1.0e-7) && (count < 5000)) {
    count++;
    std::cout << "count: " << count << std::endl;
    // update μ(y*)
    arma::colvec functand = eOfZ( designMatrixX, SMatrix, muOptimalLambda, sigmaOptimalLambda ) * muOptimalAlpha + A * muOptimalBeta;
    for (int index = 0; index < nOfData; index++) {
      if (yResponse(index) == 1.0) {
        muOptimalResponseStar(index) = functand(index) + ((1.0 / std::sqrt(2.0 * M_PI)) * std::exp(-0.5 * functand(index) * functand(index))) / (0.5 * std::erfc(-functand(index) * M_SQRT1_2));
      } else {
        muOptimalResponseStar(index) = functand(index) + ((1.0 / std::sqrt(2.0 * M_PI)) * std::exp(-0.5 * functand(index) * functand(index))) / (0.5 * std::erfc(-functand(index) * M_SQRT1_2) - 1.0);
      }
    }
    
    // update μ(α), Σ(α)
    arma::mat tempDiagonalMatrix(2 * cutOff, 2 * cutOff); tempDiagonalMatrix.eye();
    sigmaOptimalAlpha = ( m * std::exp(logH(2.0 * m - 4.0, C_sigma, std::pow(A_sigma, 2.0)) - logH(2.0 * m - 2.0, C_sigma, std::pow(A_sigma, 2.0))) * tempDiagonalMatrix + eOfZTZ( designMatrixX, SMatrix, muOptimalLambda, sigmaOptimalLambda ) ).i();
    muOptimalAlpha    = sigmaOptimalAlpha * ( eOfZ( designMatrixX, SMatrix, muOptimalLambda, sigmaOptimalLambda ).t() * muOptimalResponseStar - eOfZ( designMatrixX, SMatrix, muOptimalLambda, sigmaOptimalLambda ).t() * A * muOptimalBeta );

    // update μ(β), Σ(β)
    sigmaOptimalBeta = (priorSigmaBeta.i() + A.t() * A).i();
    muOptimalBeta = sigmaOptimalBeta * ( arma::solve( priorSigmaBeta, priorMuBeta ) + A.t() * muOptimalResponseStar - A.t() * eOfZ( designMatrixX, SMatrix, muOptimalLambda, sigmaOptimalLambda ) * muOptimalAlpha );

    // update C(σ)
    C_sigma = 0.5 * m * (arma::trace(sigmaOptimalAlpha) + arma::dot(muOptimalAlpha, muOptimalAlpha));

    // update μ(λ), Σ(λ)
    arma::mat F1(dimension, dimension); F1.zeros();
    arma::mat F2(dimension, dimension); F2.zeros();
    arma::colvec F3(dimension); F3.zeros();
    arma::colvec F4(dimension); F4.zeros();
    arma::colvec t_ij(dimension);
    arma::colvec tPlus_ijr(dimension);
    arma::colvec tMinus_ijr(dimension);
    arma::colvec AmuOptimalBeta = A * muOptimalBeta;
    arma::mat ABD = sigmaOptimalAlpha + muOptimalAlpha * muOptimalAlpha.t();
    arma::mat Amat = ABD.submat(0, 0, cutOff - 1, cutOff - 1);
    arma::mat Bmat = ABD.submat(cutOff , 0, 2 * cutOff - 1, cutOff - 1);
    arma::mat Dmat = ABD.submat(cutOff , cutOff, 2 * cutOff - 1, 2 * cutOff - 1);
    for (int i = 0; i < nOfData; i++) {
      double temp = muOptimalResponseStar(i) - AmuOptimalBeta(i);
      for (int j = 0; j < cutOff; j++) {
        t_ij = (SMatrix.row(j) % designMatrixX.row(i)).t();
        F1 += temp * std::exp( -0.5 * arma::dot(t_ij, sigmaOptimalLambda * t_ij) ) * (muOptimalAlpha(j) * cos( arma::dot(t_ij, muOptimalLambda) ) + muOptimalAlpha(j + cutOff) * sin( arma::dot( t_ij, muOptimalLambda ) )) * (t_ij * t_ij.t());
        F3 += temp * std::exp( -0.5 * arma::dot(t_ij, sigmaOptimalLambda * t_ij) ) * (muOptimalAlpha(j + cutOff) * cos( arma::dot(t_ij, muOptimalLambda) ) - muOptimalAlpha(j) * sin( arma::dot( t_ij, muOptimalLambda ) )) * t_ij;
        for (int r = 0; r < cutOff; r++) {
          tPlus_ijr = (SMatrix.row(j) % designMatrixX.row(i) + SMatrix.row(r) % designMatrixX.row(i)).t();
          tMinus_ijr = (SMatrix.row(j) % designMatrixX.row(i) - SMatrix.row(r) % designMatrixX.row(i)).t();

          F2 += (std::exp( -0.5 * arma::dot( tMinus_ijr, sigmaOptimalLambda * tMinus_ijr ) ) * ( (Amat(j, r) + Dmat(j, r) ) * cos( arma::dot( tMinus_ijr, muOptimalLambda ) ) + 2.0 * Bmat(j, r) * sin( arma::dot( tMinus_ijr, muOptimalLambda ) ) ) * (tMinus_ijr * tMinus_ijr.t())) + (std::exp( -0.5 * arma::dot( tPlus_ijr, sigmaOptimalLambda * tPlus_ijr ) ) * ( (Amat(j, r) - Dmat(j, r)) * cos( arma::dot(tPlus_ijr, muOptimalLambda) ) + 2.0 * Bmat(j, r) * sin( arma::dot(tPlus_ijr, muOptimalLambda) ) ) * (tPlus_ijr * tPlus_ijr.t()));
          F4 += (std::exp( -0.5 * arma::dot( tMinus_ijr, sigmaOptimalLambda * tMinus_ijr ) ) * ( 2.0 * Bmat(j, r) * cos( arma::dot(tMinus_ijr, muOptimalLambda) ) - (Amat(j, r) + Dmat(j, r)) * sin( arma::dot( tMinus_ijr, muOptimalLambda ) ) ) * tMinus_ijr) + std::exp( -0.5 * arma::dot( tPlus_ijr, sigmaOptimalLambda * tPlus_ijr ) ) * ( 2.0 * Bmat(j, r) * cos( arma::dot(tPlus_ijr, muOptimalLambda) ) + (Dmat(j, r) - Amat(j, r)) * sin( arma::dot( tPlus_ijr, muOptimalLambda ) ) ) * tPlus_ijr;
        }
      }
    }
    F2 *= -0.25;
    F4 *= -0.25;

    arma::mat    F5 = priorSigmaLambda.i() + F1 + F2;
    arma::colvec F6 = arma::solve(priorSigmaLambda, priorMuLambda - muOptimalLambda) + F3 + F4;
    // Adaptive nonconjugate message passing step 2
    
    double diff = -1.0;
    double tempOldLowerBound = -1.0e7;
    double tempNewLowerBound;
    double stepSize = 1.0;
    double rho = 1.5;
    double tempStepSize;
    while (diff < 0.0) {
      arma::umat isPositiveDefinite = arma::zeros<arma::umat>(dimension);
      arma::umat isSymmetric = arma::zeros<arma::umat>(dimension, dimension);
      while ((arma::accu(isSymmetric) != (dimension * dimension)) || (arma::accu(isPositiveDefinite) != dimension)) {
        count2++;
        sigmaOptimalLambda = (((1.0 - stepSize) * priorSigmaLambda.i()) + stepSize * F5).i();
        arma::cx_vec eigVal = arma::eig_gen(sigmaOptimalLambda);
        arma::colvec realEigVal = arma::real(eigVal);
        isPositiveDefinite = realEigVal > 0.0;
        isSymmetric = sigmaOptimalLambda == sigmaOptimalLambda.t();
        stepSize /= rho;
      }
      muOptimalLambda   += stepSize * sigmaOptimalLambda * F6;
      tempNewLowerBound = lowerBound(yResponse, designMatrixX, SMatrix, muOptimalLambda, sigmaOptimalLambda, muOptimalAlpha, sigmaOptimalAlpha, muOptimalBeta, sigmaOptimalBeta, A, priorMuBeta, priorSigmaBeta, priorMuLambda, priorSigmaLambda, A_sigma, C_sigma);
      diff = tempNewLowerBound - tempOldLowerBound;
      tempOldLowerBound = tempNewLowerBound;
      tempStepSize = stepSize;
      stepSize = 1.0;
      std::cout << "Count2: " << count2 << std::endl;
    }
    count2 = 0;
    stepSize = tempStepSize * rho;

    newLowerBound = lowerBound(yResponse, designMatrixX, SMatrix, muOptimalLambda, sigmaOptimalLambda, muOptimalAlpha, sigmaOptimalAlpha, muOptimalBeta, sigmaOptimalBeta, A, priorMuBeta, priorSigmaBeta, priorMuLambda, priorSigmaLambda, A_sigma, C_sigma);

    std::cout << "newLowerBound: " << newLowerBound << std::endl;
    dif = std::abs(newLowerBound - oldLowerBound);
    oldLowerBound = newLowerBound;
  }
  arma::umat yStarResult = (muOptimalResponseStar > 0.0);
  for (int i = 0; i < nOfData; i++) {
    yResponse(i) = static_cast<unsigned int>(yResponse(i));
  }
  arma::umat crossVal = yStarResult == yResponse;
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
  std::cout << "Number of iterations" << std::endl;
  std::cout << count << std::endl;
  std::cout << "Does y* belong to the first category?" << std::endl;
  std::cout << arma::vectorise(yStarResult) << std::endl;
  std::cout << "How many did this model predict?" << std::endl;
  std::cout << crossVal << std::endl;
  std::cout << arma::accu(crossVal) << std::endl;
}
