require(fAsianOptions)
require(Rcpp)
require(RcppArmadillo)
sourceCpp('auxiliary.cpp')

# Auxiliary functions
intSigma <- function(p,q,r) {
  0.5 * r^(-p/2-1)*(q * gamma(p/2+1) * Re(kummerM(q^2/(4*r), p/2+1, 3/2)) + sqrt(r)*gamma((p+1)/2) * Re(kummerM(q^2/(4*r), (p+1)/2, 1/2)))
}

E1overSigma <- function(a, b, c) {
  intSigma(2*a-1, b, c) / intSigma(2*a-2, b, c)
}

E1overSigma2 <- function(a, b, c) {
  intSigma(2*a, b, c) / intSigma(2*a-2, b, c)
}

dnormOverPnorm <- function(x) sqrt(2/pi) * exp(-x^2/2) / (erf(x/sqrt(2)) + 1)
dnormOverPnorm_Laurent <- function(x) -x - 1/x + 2/x^3 - 10/x^5
dnormOverPnormMinusOne <- function(x) exp(-x^2/2) / (sqrt(pi/2) * (erf(x/sqrt(2))+1) - 1)

Qj <- function(muPsi, sigmaPsi2, j) {
  sigmaPsi <- sqrt(sigmaPsi2)
  exp(0.5 * sigmaPsi2 * j^2 + muPsi * j) * (1 - pnorm(-muPsi/sigmaPsi - sigmaPsi * j)) + exp(0.5 * sigmaPsi2 * j^2 - muPsi * j) * (1 - pnorm(muPsi/sigmaPsi - sigmaPsi * j))
}

S1 <- function(muPsi, sigmaPsi2, w0) {
  sigmaPsi <- sqrt(sigmaPsi2)
  -w0 * (sigmaPsi * sqrt(2 / pi) * exp(-muPsi^2 / (2 * sigmaPsi2)) + muPsi * (1 - 2 * pnorm(-muPsi / sigmaPsi)))
}

S2 <- function(muPsi, sigmaPsi2, muThetaq, SigmaThetaq, a, b, c, rqt, sqt, w0, J) {
  -0.5 * E1overSigma(a, b, c) * rqt/sqt * sum((diag(SigmaThetaq) + muThetaq^2) * Qj(muPsi, sigmaPsi2, 1:J)) - 0.25 * J * (J + 1) / w0 * S1(muPsi, sigmaPsi2, w0)
}


DS1_sigmaPsi2 <- function(muPsi, sigmaPsi2, w0) {
  sigmaPsi <- sqrt(sigmaPsi2)
  -w0 * ((1/(2 * sigmaPsi) + 0.5 * sigmaPsi * muPsi^2 / sigmaPsi2^2) * sqrt(2 / pi) * exp(-0.5 * muPsi^2 / sigmaPsi2) - muPsi^2/sigmaPsi^3 * dnorm(-muPsi/sigmaPsi))
}

DS1_muPsi <- function(muPsi, sigmaPsi2, w0) {
  sigmaPsi <- sqrt(sigmaPsi2)
  -w0 * (-muPsi/sigmaPsi * sqrt(2/pi) * exp(-0.5 * muPsi^2 / sigmaPsi2) + (1 - 2 * pnorm(-muPsi/sigmaPsi)) + 2 * muPsi / sigmaPsi * dnorm(-muPsi/sigmaPsi))
}

DQj_muPsi <- function(muPsi, sigmaPsi2, j) {
  sigmaPsi <- sqrt(sigmaPsi2)
  term1 <- 0.5 * sigmaPsi2 * j^2
  term2 <- muPsi * j
  term3 <- muPsi / sigmaPsi
  term4 <- sigmaPsi * j
  term5 <- 1/sigmaPsi
  exp(term1 + term2) * dnorm(-term3 - term4) * term5 + j * exp(term1 + term2) * (1 - pnorm(-term3 - term4)) - exp(term1 - term2) * dnorm(term3 - term4) * term5 - j * exp(term1 - term2) * (1 - pnorm(term3 - term4))
}

DQj_sigmaPsi2 <- function(muPsi, sigmaPsi2, j) {
  sigmaPsi <- sqrt(sigmaPsi2)
  term1 <- 0.5 * sigmaPsi2 * j^2
  term2 <- muPsi * j
  term3 <- muPsi / sigmaPsi
  term4 <- sigmaPsi * j
  term5 <- 0.5 * muPsi / sigmaPsi^3
  term6 <- 0.5 * j / sigmaPsi
  term7 <- 0.5 * j^2
  -exp(term1 + term2) * dnorm(-term3 - term4) * (term5 - term6) + term7 * exp(term1 + term2) * (1 - pnorm(-term3 - term4)) + exp(term1 - term2) * dnorm(term3 - term4) * (term5 + term6) + term7 * exp(term1 - term2) * (1 - pnorm(term3 - term4))
}

DS2_sigmaPsi2 <- function(muPsi, sigmaPsi2, w0, a, b, c, rqt_over_sqt, SigmaThetaq, muThetaq) {
  J <- length(muThetaq) - 1
  SigmaThetaqStar <- SigmaThetaq[-1,]
  SigmaThetaqStar <- SigmaThetaqStar[,-1]
  muThetaqStar <- muThetaq[-1]
  -0.25 * J * (J + 1) / w0 * DS1_sigmaPsi2(muPsi, sigmaPsi2, w0) - 0.5 * E1overSigma(a, b, c) * rqt_over_sqt * sum((diag(SigmaThetaqStar) + muThetaqStar^2) * DQj_sigmaPsi2(muPsi, sigmaPsi2, 1:J))
}

DS2_muPsi <- function(muPsi, sigmaPsi2, w0, a, b, c, rqt_over_sqt, SigmaThetaq, muThetaq) {
  J <- length(muThetaq) - 1
  SigmaThetaqStar <- SigmaThetaq[-1,]
  SigmaThetaqStar <- SigmaThetaqStar[,-1]
  muThetaqStar <- muThetaq[-1]
  -0.25 * J * (J + 1) / w0 * DS1_muPsi(muPsi, sigmaPsi2, w0) - 0.5 * E1overSigma(a, b, c) * rqt_over_sqt * sum((diag(SigmaThetaqStar) + muThetaqStar^2) * DQj_muPsi(muPsi, sigmaPsi2, 1:J))
}

MDarrayMult <- function(MDarray, VorM) {
  ###############################################################
  ## the number of matrices inside the multi-dimensional array ##
  ## we will only be using 3-dimensional arrays #################
  ###############################################################
  d <- dim(MDarray)[3]
  temp <- array(0, dim = dim(MDarray))
  for (i in 1:d) {
    temp[,,i] <- MDarray[,,i] %*% VorM
  }
  temp
}

sumVarQuad <- function(MDarray, meanVec, covMat) {
  d <- dim(MDarray)[3]
  temp <- 0
  for (i in 1:d) {
    tempVec <- MDarray[,,i] %*% meanVec
    temp <- temp + 2 * (sum(MDarray[,,i] * covMat))^2 + 4 * sum(covMat * tcrossprod(tempVec))
  }
  temp
}

DS1_SigmaThetaq <- function(a, b, c, muPsi, sigmaPsi2, sigma0, J, rqt, sqt) {
  -0.5 * E1overSigma(a, b, c) * diag(c(1/sigma0, rqt/sqt *Qj(muPsi, sigmaPsi2, 1:J)))
}

DS2_SigmaThetaq <- function(varphi, muThetaq, SigmaThetaq, W, muYStarq, muBetaq, a, b, c, delta) {
  J <- dim(varphi)[1] - 1
  n <- dim(varphi)[3]
  temp <- matrix(0, nrow = J+1, ncol = J+1)
  for (i in 1:n) {
    tempVec <- varphi[,,i] %*% muThetaq
    temp <- temp + (6 * varphi[,,i] %*% SigmaThetaq %*% varphi[,,i] + 4 * tcrossprod(tempVec)) -2 * delta * (muYStarq[i] - sum(W[i,] * muBetaq) - delta * sum(varphi[,,i] * tcrossprod(muThetaq))) * varphi[,,i]
  }
  -0.5 * E1overSigma2(a, b, c) * temp
}

# DS2_SigmaThetaq <- function(varphi, muThetaq, SigmaThetaq, W, muYStarq, muBetaq, a, b, c, delta) {
#   J <- dim(varphi)[1] - 1
#   n <- dim(varphi)[3]
#   temp <- matrix(0, nrow = J+1, ncol = J+1)
#   for (i in 1:n) {
#     tempVec <- varphi[,,i] %*% muThetaq
#     temp <- temp + 4 * (varphi[,,i] %*% SigmaThetaq %*% varphi[,,i] + tcrossprod(tempVec)) -2 * delta * (muYStarq[i] - sum(W[i,] * muBetaq) - delta * sum(SigmaThetaq * varphi[,,i]) - delta * sum(varphi[,,i] * tcrossprod(muThetaq))) * varphi[,,i]
#   }
#   -0.5 * E1overSigma2(a, b, c) * temp
# }

DS2_muThetaq <- function(varphi, muThetaq, SigmaThetaq, W, muYStarq, muBetaq, a, b, c, delta) {
  temp = rep(0, length(muThetaq))
  n = length(muYStarq)
  for (i in 1:n) {
    tempVec = varphi[,,i] %*% muThetaq
    temp = temp + 8 * (varphi[,,i] %*% (SigmaThetaq %*% tempVec)) - 4 * delta * (muYStarq[i] - sum(W[i,] * muBetaq) - delta * sum(varphi[,,i] * SigmaThetaq) - delta * c(crossprod(muThetaq, tempVec))) * tempVec
  }
  -0.5 * E1overSigma2(a, b, c) * temp
}

DS1_muThetaq <- function(a, b, c, muPsi, sigmaPsi2, muThetaq, sigma0, rqt, sqt) {
  J <- length(muThetaq) - 1
  upsilon <- c(1/sigma0, rqt/sqt * Qj(muPsi, sigmaPsi2, 1:J))
  -E1overSigma(a, b, c) * (upsilon * muThetaq)
}

LB <- function(y, varphi, W, delta, SigmaBeta0, SigmaBetaq, muThetaq, SigmaThetaq, muYStar, a, b, c, sigma0, muPsi, sigmaPsi2, w0, muBetaq, muBeta0, r0t_half, s0t_half, r0s_half, s0s_half, rqt_half, sqt_half) {
  J <- dim(varphi)[1] - 1

  SigmaBeta0_inv <- solve(SigmaBeta0)
  upsilon <- c(1/sigma0, rqt_half/sqt_half * Qj(muPsi, sigmaPsi2, 1:J))
  
  # -0.5 * (sum(crossprod(W) * SigmaBetaq) + delta^2 * sumVarQuad(varphi, muThetaq, SigmaThetaq)) + sum(log(pnorm(muYStar)^y * (1 - pnorm(muYStar))^(1-y))) + temp  - 0.5 * ((E1overSigma(a, b, c) * sum((diag(SigmaThetaq) + muThetaq^2) * upsilon)) - determinant(SigmaThetaq)$modulus[1] - determinant(SigmaBetaq)$modulus[1] + E1overSigma2(a, b, c) * (sum(SigmaBeta0_inv * tcrossprod(muBetaq - muBeta0)) + sum(SigmaBeta0_inv * SigmaBetaq)) - log(2*pi*sigmaPsi2) + 1) + log(w0/2) + S1(muPsi, sigmaPsi2, w0) + r0s_half * log(s0s_half) - lgamma(r0s_half) - s0s_half * E1overSigma2(a, b, c) - b * E1overSigma(a, b, c) + c * E1overSigma2(a, b, c) - r0t_half * log(sqt_half) + (r0t_half - rqt_half) * digamma(rqt_half) + (1 - s0t_half/sqt_half)*rqt_half + lgamma(rqt_half)
  -0.5 * (sum(crossprod(W) * SigmaBetaq) + delta^2 * sumVarQuad(varphi, muThetaq, SigmaThetaq)) + sum(log(pnorm(muYStar)^y * (1 - pnorm(muYStar))^(1-y)))  - 0.5 * ((E1overSigma(a, b, c) * sum((diag(SigmaThetaq) + muThetaq^2) * upsilon)) - determinant(SigmaThetaq)$modulus[1] - determinant(SigmaBetaq)$modulus[1] + E1overSigma2(a, b, c) * (sum(SigmaBeta0_inv * tcrossprod(muBetaq - muBeta0)) + sum(SigmaBeta0_inv * SigmaBetaq)) - log(2*pi*sigmaPsi2) + 1) + log(w0/2) + S1(muPsi, sigmaPsi2, w0) + r0s_half * log(s0s_half) - lgamma(r0s_half) - s0s_half * E1overSigma2(a, b, c) - b * E1overSigma(a, b, c) + c * E1overSigma2(a, b, c) - r0t_half * log(sqt_half) + (r0t_half - rqt_half) * digamma(rqt_half) + (1 - s0t_half/sqt_half)*rqt_half + lgamma(rqt_half) + r0t_half * log(s0t_half) - lgamma(r0t_half)
}

# setVarPhi <- function(x, J) {
#   n <- length(x)
#   varphi <- array(0, dim = c(J+1, J+1, n))  
#   for (i in 1:n) {
#     varphi[1,1,i] <- x[i] - 0.5
#     for (j in 2:(J+1)) {
#       varphi[1,j,i] <- varphi[j,1,i] <- sqrt(2)/(pi * j) * sin(pi * j * x[i]) - sqrt(2)/((pi * j)^2) * (1 - cos(pi * j))
#       for (k in j:(J+1)) {
#         if (k == j) {
#           varphi[j,j,i] <- sin(2*pi*j*x[i])/(2*pi*j) + x[i] - 0.5
#         } else {
#           varphi[j,k,i] <- varphi[k,j,i] <- sin(pi*(j+k)*x[i])/(pi*(j+k)) + sin(pi*(j-k)*x[i])/(pi*(j-k)) - (1-cos(pi*(j+k)))/((pi*(j+k))^2) - (1-cos(pi*(j-k)))/((pi*(j-k))^2)
#         }
#       }
#     }
#   }
#   varphi
# }

# setVarPhi2 <- function(x.grid, J) {
#   n.grid <- length(x.grid)
#   dmatsfull<-array(0,dim=c(n.grid,J+1,J+1))
#   dmatsfull[,1,1]<-x.grid-0.5
#   temparg1<-outer(x.grid,(1:J))
#   temparg2<-matrix(rep((1:J),each=n.grid),nrow=n.grid,byrow=FALSE)
#   dmatsfull[,1,(2:(J+1))]<-sqrt(2)/(pi*temparg2)*sin(pi*temparg1)-sqrt(2)/((pi*temparg2)^2)*(1-cos(pi*temparg2))
#   dmatsfull[,(2:(J+1)),1]<-dmatsfull[,1,(2:(J+1))]
#   temparg1<-outer((1:J),(1:J),'+')
#   temparg2<-outer(rep(1,times=n.grid),temparg1)
#   temparg1<-outer(x.grid,temparg1)
#   dmatsfull[,(2:(J+1)),(2:(J+1))]<-sin(pi*temparg1)/(pi*temparg2)-(1-cos(pi*temparg2))/(pi*temparg2)^2
#   temparg1<-outer((1:J),(1:J),'-')
#   temparg2<-outer(rep(1,times=n.grid),temparg1)
#   temparg1<-outer(x.grid,temparg1)
#   dmatsfull[,(2:(J+1)),(2:(J+1))]<-dmatsfull[,(2:(J+1)),(2:(J+1))]+sin(pi*temparg1)/(pi*temparg2)-(1-cos(pi*temparg2))/(pi*temparg2)^2
#   temparg1<-outer(x.grid,(1:J))
#   temparg2<-outer(rep(1,times=n.grid),(1:J))
#   temparg3<-outer(x.grid,rep(1,times=J))
#   tempmat<-sin(2*pi*temparg1)/(2*pi*temparg2)+temparg3-0.5
#   for(j in 2:(J+1)) {
#     dmatsfull[,j,j]<-tempmat[,(j-1)]
#   }
#   dmatsfull
# }


VB <- function(y, x, W, J, delta, r0s, s0s, r0t, s0t, w0, muBeta0, SigmaBeta0, muTheta0, SigmaTheta0, sigma0, tol = 1.0e-04) {
  ## Set design matrix varphi
  n <- length(y)
  p <- dim(W)[2]
  varphi <- setVarPhi(x, J)
  if(!is.matrix(W)) as.matrix(W)
  WtW <- crossprod(W)
  SigmaTheta0_inv <- solve(SigmaTheta0)
  SigmaBeta0_inv <- solve(SigmaBeta0)
  SigmaBeta0_inv_muBeta0 <- SigmaBeta0_inv %*% muBeta0

  r0s_half <- 0.5 * r0s
  s0s_half <- 0.5 * s0s
  r0t_half <- 0.5 * r0t
  s0t_half <- 0.5 * s0t

  ## Initialize variational parameters
  rqt <- r0t + J
  sqt <- s0t + 10
  rqt_half <- 0.5 * rqt
  sqt_half <- 0.5 * sqt
  rqt_over_sqt <- rqt/sqt
  muBetaq <- rep(0.1, p)
  SigmaBetaq <- diag(0.1, p)
  muThetaq <- rep(0.1, J+1)
  SigmaThetaq <- diag(0.1, J+1)
  SigmaThetaOldInv <- diag(J+1)
  muPsi <- 1
  sigmaPsi2 <- 1/(J*J)
  sigmaPsi <- sqrt(sigmaPsi2)
  muYStar <- rep(0.1, n)
  muYStarq <- rep(0.1, n)
  a <- 0.25 * (J+1) + r0s_half + 0.5 * p + 1
  b <- 1
  c <- 10
  lbold <- LB(y, varphi, W, delta, SigmaBeta0, SigmaBetaq, muThetaq, SigmaThetaq, muYStar, a, b, c, sigma0, muPsi, sigmaPsi2, w0, muBetaq, muBeta0, r0t_half, s0t_half, r0s_half, s0s_half, rqt_half, sqt_half)
  cat('first LB: ', lbold, '\n')
  lbnew <- 0
  dif <- tol + 1
  lbrecord <- c(lbold)
  count <- 0
  while(dif > tol) {
    count <- count + 1
    
    # Update theta

    ## SigmaTheta
    step <- 1
    SigmaThetaNewInv <- -2 * DS1_SigmaThetaq(a, b, c, muPsi, sigmaPsi2, sigma0, J, rqt, sqt) + DS2_SigmaThetaq(varphi, muThetaq, SigmaThetaq, W, muYStarq, muBetaq, a, b, c, delta)
    SigmaThetaOldInvTry <- (1-step) * SigmaThetaOldInv + step * SigmaThetaNewInv
    SigmaThetaTry <- solve(SigmaThetaOldInvTry)
    EV <- eigen(SigmaThetaTry, only.values = TRUE)$values

    while(isSymmetric(SigmaThetaTry, tol = 1.0e-10) == FALSE | any(Re(EV) < 0)) {
      step <- 2/3 * step
      SigmaThetaOldInvTry <- (1-step) * SigmaThetaOldInv + step * SigmaThetaNewInv
      SigmaThetaTry <- solve(SigmaThetaOldInvTry)
      EV <- eigen(SigmaThetaTry, only.values = TRUE)$values
    }
    SigmaThetaOldInv <- SigmaThetaOldInvTry
    SigmaThetaq <- SigmaThetaTry

    ## muTheta
    muThetaq <- muThetaq + step * SigmaThetaq %*% (DS1_muThetaq(a, b, c, muPsi, sigmaPsi2, muThetaq, sigma0, rqt, sqt) + DS2_muThetaq(varphi, muThetaq, SigmaThetaq, W, muYStarq, muBetaq, a, b, c, delta))
    cat("DS1_muThetaq: ", DS1_muThetaq(a, b, c, muPsi, sigmaPsi2, muThetaq, sigma0, rqt, sqt), '\n')
    cat("DS2_muThetaq: ", DS2_muThetaq(varphi, muThetaq, SigmaThetaq, W, muYStarq, muBetaq, a, b, c, delta), '\n')
    # SigmaThetaqNew <- solve((1-step) * solve(SigmaThetaq) + step * temp1)
    # EV <- eigen(SigmaThetaqNew, only.values = TRUE)$values

    # while(isSymmetric(SigmaThetaqNew, tol=1.0e-10)==FALSE | any(Re(EV) < 0)) {
    #   step <- 2/3 * step
    #   SigmaThetaqNew <- solve((1-step) * solve(SigmaThetaq) + step * temp1)
    #   EV <- eigen(SigmaThetaqNew, only.values=TRUE)$values
    # }

    # SigmaThetaq <- SigmaThetaqNew
    # muThetaq <- muThetaq + step*SigmaThetaq %*% (DS1_muThetaq(a, b, c, muPsi, sigmaPsi2, muThetaq, sigma0, rqt, sqt) + DS2_muThetaq(varphi, muThetaq, SigmaThetaq, W, muYStarq, muBetaq, a, b, c, delta))
    
    # Update psi
    step_psi = 1
    sigmaPsi2_old = sigmaPsi2
    muPsi_old = muPsi
    sigmaPsi2_new = -0.5 / (DS1_sigmaPsi2(muPsi, sigmaPsi2, w0) + DS2_sigmaPsi2(muPsi, sigmaPsi2, w0, a, b, c, rqt_over_sqt, SigmaThetaq, muThetaq))
    temp2 = sigmaPsi2_new
    muPsi_old = muPsi

    temp = DS1_muPsi(muPsi, sigmaPsi2, w0) + DS2_muPsi(muPsi, sigmaPsi2, w0, a, b, c, rqt_over_sqt, SigmaThetaq, muThetaq)
    while (sigmaPsi2_new < 0) {
      step_psi <- 2/3 * step_psi
      sigmaPsi2_new <- 1 / (1/sigmaPsi2_old + step_psi * (1/temp2 - 1/sigmaPsi2_old))
    }
    sigmaPsi2 <- sigmaPsi2_new
    # while (sigmaPsi2_new < 0)
    # {
    #   step_psi = 2/3 * step_psi
    #   sigmaPsi2_new = 1 / (1/sigmaPsi2_old + step_psi * (1/sigmaPsi2_new - 1/sigmaPsi2_old))
    # }
    # sigmaPsi2 = sigmaPsi2_new
    sigmaPsi = sqrt(sigmaPsi2)
    muPsi = muPsi + step_psi * sigmaPsi2 * temp

    # cat("Trace of SigmaBetaq: ", sum(diag(SigmaBetaq)), '\n')
    # cat("Sum of squared muBetaq: ", sum(muBetaq^2), '\n')
    # cat("E1overSigma2: ", E1overSigma2(a, b, c), '\n')
    # cat("Wt %*% muYStarq: ", c(crossprod(W, muYStarq)), '\n')

    # Update tau^2
    SigmaThetaqStar <- SigmaThetaq[-1,]
    SigmaThetaqStar <- SigmaThetaqStar[,-1]
    muThetaqStar <- muThetaq[-1]
    sqt <- s0t + E1overSigma(a, b, c) * sum((diag(SigmaThetaqStar) + muThetaqStar^2) * Qj(muPsi, sigmaPsi2, 1:J))
    sqt_half <- 0.5 * sqt
    rqt_over_sqt <- rqt/sqt
    cat("sqt: ", sqt, '\n')
    cat("rqt_over_sqt: ", rqt_over_sqt, '\n')
    cat("rqt_over_sqt * Qj: ", rqt_over_sqt * Qj(muPsi, sigmaPsi2, 1:J), '\n')
    cat("E1overSigma: ", E1overSigma(a, b, c), '\n')
    # Update sigma^2
    b <- -0.5 * sum((diag(SigmaThetaq) + muThetaq^2) * c(1/sigma0, rqt_over_sqt * Qj(muPsi, sigmaPsi2, 1:J)))
    c <- 0.5 * (s0s + sum(SigmaBeta0_inv * tcrossprod(muBetaq - muBeta0)) + sum(SigmaBeta0_inv * SigmaBetaq))
    cat("b: ", b, ", c: ", c, '\n')
    # cat('muBetaq: ', muBetaq, '\n')
    # cat('sum(varphi[,,1] * SigmaThetaq): ', sum(varphi[,,1] * SigmaThetaq), '\n')
    # cat('sum(varphi[,,1] * tcrossprod(muThetaq)): ', sum(varphi[,,1] * tcrossprod(muThetaq)), '\n')
    # Update y*
    for (i in 1:n) {
      tempMuYStar <- sum(W[i,] * muBetaq) + delta * sum(varphi[,,i] * SigmaThetaq) + delta * sum(varphi[,,i] * tcrossprod(muThetaq))
      if (y[i] == 1 & tempMuYStar < -8) {
        muYStar[i] <- tempMuYStar
        muYStarq[i] <- tempMuYStar + dnormOverPnorm_Laurent(tempMuYStar)
      } else if (y[i] == 1 & tempMuYStar >= -8) {
        muYStar[i] <- tempMuYStar
        muYStarq[i] <- tempMuYStar + dnormOverPnorm(tempMuYStar)
      } else {
        muYStar[i] <- tempMuYStar
        muYStarq[i] <- tempMuYStar + dnormOverPnormMinusOne(tempMuYStar)
      }
      # muYStar[i] <- sum(W[i,] * muBetaq) + delta * sum(varphi[,,i] * SigmaThetaq) + delta * sum(varphi[,,i] * tcrossprod(muThetaq))
      # cat('delta * sum(varphi[,,i] * SigmaThetaq): ', delta * sum(varphi[,,i] * SigmaThetaq), '\n')
      # cat('delta * sum(varphi[,,i] * tcrossprod(muThetaq): ', delta * sum(varphi[,,i] * tcrossprod(muThetaq)), '\n')
    }
    # plot(density(muYStar))
    # muYStarq <- muYStar + dnorm(muYStar)/((pnorm(muYStar)^y)*(pnorm(muYStar)-1)^(1-y))
    cat('sum of muYStarq: ', sum(muYStarq), '\n')
    # stop('stopped')
    # Update beta
    SigmaBetaq <- solve(E1overSigma2(a, b, c) * SigmaBeta0_inv + WtW)

    tempBeta <- rep(0, p)
    for (i in 1:n) {
      tempBeta <- tempBeta + (sum(varphi[,,i] * SigmaThetaq) + sum(varphi[,,i] * tcrossprod(muThetaq))) * W[i]
    }
    muBetaq <- SigmaBetaq %*% (E1overSigma2(a, b, c) * SigmaBeta0_inv_muBeta0 + c(crossprod(W, muYStarq)) - delta * tempBeta)

    lbnew <- LB(y, varphi, W, delta, SigmaBeta0, SigmaBetaq, muThetaq, SigmaThetaq, muYStar, a, b, c, sigma0, muPsi, sigmaPsi2, w0, muBetaq, muBeta0, r0t_half, s0t_half, r0s_half, s0s_half, rqt_half, sqt_half)
    diff <- lbnew - lbold
    cat('lbnew: ', lbnew, '\n')
    cat('lbold: ', lbold, '\n')
    lboldtry <- lbold
    if (diff < 0) {
      # Update theta
      ## SigmaTheta
      step <- 1
      SigmaThetaNewInv <- -2 * DS1_SigmaThetaq(a, b, c, muPsi, sigmaPsi2, sigma0, J, rqt, sqt) + DS2_SigmaThetaq(varphi, muThetaq, SigmaThetaq, W, muYStarq, muBetaq, a, b, c, delta)
      SigmaThetaOldInvTry <- (1-step) * SigmaThetaOldInv + step * SigmaThetaNewInv
      SigmaThetaTry <- solve(SigmaThetaOldInvTry)
      EV <- eigen(SigmaThetaTry, only.values = TRUE)$values

      while(isSymmetric(SigmaThetaTry, tol = 1.0e-10) == FALSE | any(Re(EV) < 0)) {
        step <- 2/3 * step
        SigmaThetaOldInvTry <- (1-step) * SigmaThetaOldInv + step * SigmaThetaNewInv
        SigmaThetaTry <- solve(SigmaThetaOldInvTry)
        EV <- eigen(SigmaThetaTry, only.values = TRUE)$values
      }
      SigmaThetaOldInv <- SigmaThetaOldInvTry
      SigmaThetaq <- SigmaThetaTry

      ## muTheta
      muThetaq <- muThetaq + step * SigmaThetaq %*% (DS1_muThetaq(a, b, c, muPsi, sigmaPsi2, muThetaq, sigma0, rqt, sqt) + DS2_muThetaq(varphi, muThetaq, SigmaThetaq, W, muYStarq, muBetaq, a, b, c, delta))
      cat("SigmaThetaq: ", SigmaThetaq, '\n')
      # Update psi
      step_psi = 1
      sigmaPsi2_old = sigmaPsi2
      muPsi_old = muPsi
      sigmaPsi2_new = -0.5 / (DS1_sigmaPsi2(muPsi, sigmaPsi2, w0) + DS2_sigmaPsi2(muPsi, sigmaPsi2, w0, a, b, c, rqt_over_sqt, SigmaThetaq, muThetaq))
      temp2 = sigmaPsi2_new
      muPsi_old = muPsi

      temp = DS1_muPsi(muPsi, sigmaPsi2, w0) + DS2_muPsi(muPsi, sigmaPsi2, w0, a, b, c, rqt_over_sqt, SigmaThetaq, muThetaq)
      while (sigmaPsi2_new < 0)
      {
        step_psi = 2/3 * step_psi
        sigmaPsi2_new = 1 / (1/sigmaPsi2_old + step_psi * (1/temp2 - 1/sigmaPsi2_old))
      }
      sigmaPsi2 = sigmaPsi2_new
      sigmaPsi = sqrt(sigmaPsi2)
      muPsi = muPsi + step_psi * sigmaPsi2 * temp

      # Update y*
      for (i in 1:n) {
        tempMuYStar <- sum(W[i,] * muBetaq) + delta * sum(varphi[,,i] * SigmaThetaq) + delta * sum(varphi[,,i] * tcrossprod(muThetaq))
        if (y[i] == 1 & tempMuYStar < -8) {
          muYStar[i] <- tempMuYStar
          muYStarq[i] <- tempMuYStar + dnormOverPnorm_Laurent(tempMuYStar)
        } else if (y[i] == 1 & tempMuYStar >= -8) {
          muYStar[i] <- tempMuYStar
          muYStarq[i] <- tempMuYStar + dnormOverPnorm(tempMuYStar)
        } else {
          muYStar[i] <- tempMuYStar
          muYStarq[i] <- tempMuYStar + dnormOverPnormMinusOne(tempMuYStar)
        }
        # muYStar[i] <- sum(W[i,] * muBetaq) + delta * sum(varphi[,,i] * SigmaThetaq) + delta * sum(varphi[,,i] * tcrossprod(muThetaq))
        # cat('delta * sum(varphi[,,i] * SigmaThetaq): ', delta * sum(varphi[,,i] * SigmaThetaq), '\n')
        # cat('delta * sum(varphi[,,i] * tcrossprod(muThetaq): ', delta * sum(varphi[,,i] * tcrossprod(muThetaq)), '\n')
      }
      # muYStarq <- muYStar + dnorm(muYStar)/((pnorm(muYStar)^y)*(pnorm(muYStar)-1)^(1-y))
      cat('sum of muYStarq: ', sum(muYStarq), '\n')

      # Update beta
      SigmaBetaq <- solve(E1overSigma2(a, b, c) * SigmaBeta0_inv + WtW)

      tempBeta <- rep(0, p)
      for (i in 1:n) {
        tempBeta <- tempBeta + (sum(varphi[,,i] * SigmaThetaq) + sum(varphi[,,i] * tcrossprod(muThetaq))) * W[i]
      }
      muBetaq <- SigmaBetaq %*% (E1overSigma2(a, b, c) * SigmaBeta0_inv_muBeta0 + c(crossprod(W, muYStarq)) - delta * tempBeta)
      # cat("Trace of SigmaBetaq: ", sum(diag(SigmaBetaq)), '\n')
      # cat("Sum of squared muBetaq: ", sum(muBetaq^2), '\n')
      # cat("E1overSigma2: ", E1overSigma2(a, b, c), '\n')
      # cat("Wt %*% muYStarq: ", c(crossprod(W, muYStarq)), '\n')

      # Update tau^2
      SigmaThetaqStar <- SigmaThetaq[-1,]
      SigmaThetaqStar <- SigmaThetaqStar[,-1]
      muThetaqStar <- muThetaq[-1]
      sqt <- s0t + E1overSigma(a, b, c) * sum((diag(SigmaThetaqStar) + muThetaqStar^2) * Qj(muPsi, sigmaPsi2, 1:J))
      sqt_half <- 0.5 * sqt
      rqt_over_sqt <- rqt/sqt

      # Update sigma^2
      b <- -0.5 * sum((diag(SigmaThetaq) + muThetaq^2) * c(1/sigma0, rqt_over_sqt * Qj(muPsi, sigmaPsi2, 1:J)))
      c <- 0.5 * (s0s + sum(SigmaBeta0_inv * tcrossprod(muBetaq - muBeta0)) + sum(SigmaBeta0_inv * SigmaBetaq))
      cat('b: ', b, '\n')
      cat('c: ', c, '\n')
      lbnew <- LB(y, varphi, W, delta, SigmaBeta0, SigmaBetaq, muThetaq, SigmaThetaq, muYStar, a, b, c, sigma0, muPsi, sigmaPsi2, w0, muBetaq, muBeta0, r0t_half, s0t_half, r0s_half, s0s_half, rqt_half, sqt_half)
      diff <- lbnew - lboldtry
      lbold <- lboldtry
      lboldtry <- lbnew
      # lbold <- lbnew
      cat('lbnew: ', lbnew, '\n')
      cat('lbold: ', lbold, '\n')
    }
    dif <- (lbnew - lbold)/abs(lbnew)
    lbold <- lbnew
    lbrecord <- c(lbrecord, lbnew)
    cat("count: ", count, ", lbnew: ", lbnew, ", dif: ", dif, '\n')
  }
  list(count = count, muThetaq = muThetaq, SigmaThetaq = SigmaThetaq, muBetaq = muBetaq, SigmaBetaq = SigmaBetaq, rqt = rqt, sqt = sqt, a = a, b = b, c = c, muYStarq = muYStarq, muYStar = muYStar, lbrecord = lbrecord)
}