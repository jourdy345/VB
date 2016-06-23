require(fAsianOptions)
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

Qj <- function(muPsi, sigmaPsi2, j) {
  sigmaPsi <- sqrt(sigmaPsi2)
  exp(0.5 * sigmaPsi2 * j^2 + muPsi * j) * (1 - pnorm(-muPsi/sigmaPsi - sigmaPsi * j)) + exp(9.5 * sigmaPsi2 * j^2 - muPsi * j) * (1 - pnorm(muPsi/sigmaPsi - sigmaPsi * j))
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
  term6 <- 0.5 * j / sigmaPsi2
  term7 <- 0.5 * j^2
  -exp(term1 + term2) * dnorm(-term3 - term4) * (term5 - term6) + term7 * exp(term1 + term2) * (1 - pnorm(-term3 - term4)) + exp(term1 - term2) * dnorm(term3 - term4) * (term5 + term6) + term7 * exp(term1 - term2) * (1 - pnorm(term3 - term4))
}

DS2_sigmaPsi2 <- function(muPsi, sigmaPsi2, w0, a, b, c, rqt_over_sqt, SigmaThetaq, muThetaq) {
  J <- length(muThetaq)
  -0.25 * J * (J + 1) / w0 * DS1_sigmaPsi2(muPsi, sigmaPsi2, w0) - 0.5 * E1overSigma(a, b, c) * rqt_over_sqt * sum((diag(SigmaThetaq) + muThetaq^2) * DQj_sigmaPsi2(muPsi, sigmaPsi2, 1:J))
}

DS2_muPsi <- function(muPsi, sigmaPsi2, w0, a, b, c, rqt_over_sqt, SigmaThetaq, muThetaq) {
  J <- length(muThetaq)
  -0.25 * J * (J + 1) / w0 * DS1_muPsi(muPsi, sigmaPsi2, w0) - 0.5 * E1overSigma(a, b, c) * rqt_over_sqt * sum((diag(SigmaThetaq) + muThetaq^2) * DQj_muPsi(muPsi, sigmaPsi2, 1:J))
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
    temp <- temp + 4 * (varphi[,,i] %*% SigmaThetaq %*% varphi[,,i] + tcrossprod(tempVec)) -2 * delta * (muYStarq[i] - sum(W[i,] * muBetaq) - delta * sum(SigmaThetaq * varphi[,,i]) - delta * sum(SigmaThetaq * tcrossprod(muThetaq)))
  }
  -0.5 * E1overSigma2(a, b, c) * temp
}

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
  sum(crossprod(W) * SigmaBetaq) + delta^2 * sumVarQuad(varphi, muThetaq, SigmaThetaq) + sum(log(pnorm(muYStar)^y * (1 - pnorm(muYStar))^(1-y))) - 0.5 * (E1overSigma(a, b, c) * sum((diag(SigmaThetaq) + c(crossprod(muThetaq))) * upsilon) - determinant(SigmaThetaq)$modulus[1] - determinant(SigmaBetaq)$modulus[1] + E1overSigma2(a, b, c) * (sum(SigmaBeta0_inv * tcrossprod(muBetaq - muBeta0)) + sum(SigmaBeta0_inv * SigmaBetaq)) - log(2*pi*sigmaPsi2) + 1) + log(w0/2) + S1(muPsi, sigmaPsi2, w0) + r0s_half * log(s0s_half) - lgamma(r0s_half) - s0s_half * E1overSigma2(a, b, c) - b * E1overSigma(a, b, c) + c * E1overSigma2(a, b, c) - r0t_half * log(sqt_half) + (r0t_half - rqt_half) * digamma(rqt_half) + (1 - s0t_half/sqt_half)*rqt_half + lgamma(rqt_half)
}

setVarPhi <- function(x, J) {
  n <- length(x)
  varphi = array(0, dim = c(J+1, J+1, n))  
  for (i in 1:n) {
    varphi[1,1,i] <- x[i] - 0.5
    for (j in 1:J+1) {
      varphi[1,j,i] <- varphi[j,1,i] <- sqrt(2)/(pi * j) * sin(pi * j * x[i]) - sqrt(2)/((pi * j)^2) * (1 - cos(pi * j))
      for (k in 2:J+1) {
        if (k == j) {
          varphi[j,j,i] <- sin(2*pi*j*x[i])/(2*pi*j) + x[i] - 0.5
        } else {
          varphi[j,k,i] <- sin(pi*(j+k)*x[i])/(pi*(j+k)) + sin(pi*(j-k)*x[i])/(pi*(j-k)) - (1-cos(pi*(j+k)))/((pi*(j+k))^2) - (1-cos(pi*(j-k)))/((pi*(j-k))^2)
        }
      }
    }
  }
  varphi
}

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
  sqt <- s0t
  rqt_half <- 0.5 * rqt
  sqt_half <- 0.5 * sqt
  rqt_over_sqt <- rqt/sqt
  muBetaq <- rep(0, p)
  SigmaBetaq <- diag(1, p)
  muThetaq <- rep(0, J+1)
  SigmaThetaq <- diag(1, J+1)
  muPsi <- 0.1
  sigmaPsi2 <- muPsi^2/100
  sigmaPsi <- sqrt(sigmaPsi2)
  muYStar <- rep(0, n)
  muYStarq <- rep(0, n)
  a <- 0.25 * (J+1) + r0s_half + 0.5 * p + 1
  b <- 2
  c <- 2
  lbold <- LB(y, varphi, W, delta, SigmaBeta0, SigmaBetaq, muThetaq, SigmaThetaq, muYStar, a, b, c, sigma0, muPsi, sigmaPsi2, w0, muBetaq, muBeta0, r0t_half, s0t_half, r0s_half, s0s_half, rqt_half, sqt_half)
  lbnew <- 0
  dif <- tol + 1
  lbrecord <- c(lbold)
  count <- 0
  while(dif > tol) {
    count <- count + 1
    
    # Update theta
    step <- 1
    temp1 <- SigmaTheta0_inv + DS1_SigmaThetaq(a, b, c, muPsi, sigmaPsi2, sigma0, J, rqt, sqt) + DS2_SigmaThetaq(varphi, muThetaq, SigmaThetaq, W, muYStarq, muBetaq, a, b, c, delta)
    SigmaThetaqNew <- solve((1-step) * solve(SigmaThetaq) + step * temp1)
    EV <- eigen(SigmaThetaqNew, only.values = TRUE)$values

    while(isSymmetric(SigmaThetaqNew, tol=1.0e-10)==FALSE | any(Re(EV) < 0)) {
      step <- 2/3 * step
      SigmaThetaqNew <- solve((1-step) * solve(SigmaThetaq) + step * temp1)
      EV <- eigen(SigmaThetaqNew, only.values=TRUE)$values
    }

    SigmaThetaq <- SigmaThetaqNew
    muThetaq <- muThetaq + step*SigmaThetaq %*% (DS1_muThetaq(a, b, c, muPsi, sigmaPsi2, muThetaq, sigma0, rqt, sqt) + DS2_muThetaq(varphi, muThetaq, SigmaThetaq, W, muYStarq, muBetaq, a, b, c, delta))
    
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
      sigmaPsi2_new = 1 / (1/sigmaPsi2_old + step_psi * (1/sigmaPsi2_new - 1/sigmaPsi2_old))
    }
    sigmaPsi2 = sigmaPsi2_new
    sigmaPsi = sqrt(sigmaPsi2)
    muPsi = muPsi + step_psi * sigmaPsi2 * temp

    # Update y*
    for (i in 1:n) {
      muYStar[i] <- sum(W[i,] * muBetaq) + delta * sum(varphi[,,i] * SigmaThetaq) + delta * sum(varphi[,,i] * tcrossprod(muThetaq))
    }
    muYStarq <- muYStar + dnorm(muYStar)/((pnorm(muYStar)^y)*(pnorm(muYStar)-1)^(1-y))

    # Update beta
    SigmaBetaq <- solve(E1overSigma2(a, b, c) * SigmaBeta0_inv + WtW)

    tempBeta <- rep(0, n)
    for (i in 1:n) {
      tempBeta <- tempBeta + (sum(varphi[,,i] * SigmaThetaq) + sum(varphi[,,i] * tcrossprod(muThetaq))) * W[i]
    }
    muBetaq <- SigmaBetaq %*% (E1overSigma2(a, b, c) * SigmaBeta0_inv_muBeta0 + crossprod(W, muYStarq - delta * tempBeta))

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
    cat("b: ", b, ", c: ", c, '\n')

    lbnew <- LB(y, varphi, W, delta, SigmaBeta0, SigmaBetaq, muThetaq, SigmaThetaq, muYStar, a, b, c, sigma0, muPsi, sigmaPsi2, w0, muBetaq, muBeta0, r0t_half, s0t_half, r0s_half, s0s_half, rqt_half, sqt_half)
    diff <- lbnew - lbold
    if (diff < 0) {
      # Update theta
      step <- 1
      temp1 <- SigmaTheta0_inv + DS1_SigmaThetaq(a, b, c, muPsi, sigmaPsi2, sigma0, J, rqt, sqt) + DS2_SigmaThetaq(varphi, muThetaq, SigmaThetaq, W, muYStarq, muBetaq, a, b, c, delta)
      SigmaThetaqNew <- solve((1-step) * solve(SigmaThetaq) + step * temp1)
      EV <- eigen(SigmaThetaqNew, only.values = TRUE)$values

      while(isSymmetric(SigmaThetaqNew, tol=1.0e-10)==FALSE | any(Re(EV) < 0)) {
        step <- 2/3 * step
        SigmaThetaqNew <- solve((1-step) * solve(SigmaThetaq) + step * temp1)
        EV <- eigen(SigmaThetaqNew, only.values=TRUE)$values
      }

      SigmaThetaq <- SigmaThetaqNew
      muThetaq <- muThetaq + step*SigmaThetaq %*% (DS1_muThetaq(a, b, c, muPsi, sigmaPsi2, muThetaq, sigma0, rqt, sqt) + DS2_muThetaq(varphi, muThetaq, SigmaThetaq, W, muYStarq, muBetaq, a, b, c, delta))

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
        sigmaPsi2_new = 1 / (1/sigmaPsi2_old + step_psi * (1/sigmaPsi2_new - 1/sigmaPsi2_old))
      }
      sigmaPsi2 = sigmaPsi2_new
      sigmaPsi = sqrt(sigmaPsi2)
      muPsi = muPsi + step_psi * sigmaPsi2 * temp

      # Update y*
      for (i in 1:n) {
        muYStar[i] <- sum(W[i,] * muBetaq) + delta * sum(varphi[,,i] * SigmaThetaq) + delta * sum(varphi[,,i] * tcrossprod(muThetaq))
      }
      muYStarq <- muYStar + dnorm(muYStar)/((pnorm(muYStar)^y)*(pnorm(muYStar)-1)^(1-y))

      # Update beta
      SigmaBetaq <- solve(E1overSigma2(a, b, c) * SigmaBeta0_inv + WtW)

      tempBeta <- rep(0, n)
      for (i in 1:n) {
        tempBeta <- tempBeta + (sum(varphi[,,i] * SigmaThetaq) + sum(varphi[,,i] * tcrossprod(muThetaq))) * W[i]
      }
      muBetaq <- SigmaBetaq %*% (E1overSigma2(a, b, c) * SigmaBeta0_inv_muBeta0 + crossprod(W, muYStarq - delta * tempBeta))

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

      lbnew <- LB(y, varphi, W, delta, SigmaBeta0, SigmaBetaq, muThetaq, SigmaThetaq, muYStar, a, b, c, sigma0, muPsi, sigmaPsi2, w0, muBetaq, muBeta0, r0t_half, s0t_half, r0s_half, s0s_half, rqt_half, sqt_half)
      diff <- lbnew - lbnew
    }
    dif <- (lbnew - lbold)/abs(lbnew)
    lbold = lbnew
    lbrecord <- c(lbrecord, lbnew)
    cat("count: ", count, ", lbnew: ", lbnew, ", dif: ", dif, '\n')
  }
  list(muThetaq = muThetaq, SigmaThetaq = SigmaThetaq, muBetaq = muBetaq, SigmaBetaq = SigmaBetaq, rqt = rqt, sqt = sqt, a = a, b = b, c = c, muYStarq = muYStarq, muYStar = muYStar, lbrecord = lbrecord)
}