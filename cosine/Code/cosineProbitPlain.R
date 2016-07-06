# require(Rcpp)
# require(RcppArmadillo)
# sourceCpp('auxiliary.cpp')

# dnormOverPnorm <- function(x) sqrt(2 / pi) * exp(-x^2 / 2) / (erf(x / sqrt(2)) + 1)
# dnormOverPnorm_Laurent <- function(x) -x - 1 / x + 2 / x^3 - 10 / x^5
# dnormOverPnormMinusOne <- function(x) exp(-x^2 / 2) / (sqrt(pi / 2) * (erf(x / sqrt(2)) + 1) - 1)

Qj <- function(muPsi, sigmaPsi2, sigmaPsi, j) {
  t1 <- 0.5 * sigmaPsi2 * j^2
  t2 <- muPsi * j
  t3 <- muPsi / sigmaPsi
  t4 <- sigmaPsi * j
  exp(t1 + t2) * (1 - pnorm(-t3 - t4)) + exp(t1 - t2) * (1 - pnorm(t3 - t4))
}

S1 <- function(muPsi, sigmaPsi2, sigmaPsi, w0) {
  -w0 * (sigmaPsi * sqrt(2/pi) * exp(-0.5 * muPsi^2 / sigmaPsi2) + muPsi * (1 - 2 * pnorm(-muPsi / sigmaPsi)))
}

S2 <- function(muPsi, sigmaPsi2, sigmaPsi, w0, rqt, sqt, SigmaThetaq, muThetaq, J) {
  -0.5 * rqt/sqt * sum((diag(SigmaThetaq) + muThetaq^2) * Qj(muPsi, sigmaPsi2, sigmaPsi, 1:J)) - 0.25 * J * (J + 1) / w0 * S1(muPsi, sigmaPsi2, sigmaPsi, w0)
}

DQj_muPsi <- function(muPsi, sigmaPsi2, sigmaPsi, j) {
  t1 <- 0.5 * sigmaPsi2 * j^2
  t2 <- muPsi * j
  t3 <- muPsi / sigmaPsi
  t4 <- sigmaPsi * j
  t5 <- 1 / sigmaPsi
  exp(t1 + t2) * dnorm(-t3 - t4) * t5 + j * exp(t1 + t2) * (1 - pnorm(-t3 - t4)) - exp(t1 - t2) * dnorm(t3 - t4) * t5 - j * exp(t1 - t2) * (1 - pnorm(t3 - t4))
}

DQj_sigmaPsi2 <- function(muPsi, sigmaPsi2, sigmaPsi, j) {
  t1 <- 0.5 * sigmaPsi2 * j^2
  t2 <- muPsi * j
  t3 <- muPsi / sigmaPsi
  t4 <- sigmaPsi * j
  t5 <- 0.5 * muPsi / sigmaPsi^3
  t6 <- 0.5 * j / sigmaPsi
  t7 <- 0.5 * j^2
  -exp(t1 + t2) * dnorm(-t3 - t4) * (t5 - t6) + t7 * exp(t1 + t2) * (1 - pnorm(-t3 - t4)) + exp(t1 - t2) * dnorm(t3 - t4) * (t5 + t6) + t7 * exp(t1 - t2) * (1 - pnorm(t3 - t4))
}

DS1_muPsi <- function(muPsi, sigmaPsi2, sigmaPsi, w0) {
  -w0 * ((0.5 / sigmaPsi * (1 + muPsi^2 / sigmaPsi2)) * sqrt(2 / pi) * exp(-0.5 * muPsi^2 / sigmaPsi2) - muPsi^2 / sigmaPsi^3 * dnorm(-muPsi / sigmaPsi))
}

DS1_sigmaPsi2 <- function(muPsi, sigmaPsi2, sigmaPsi, w0) {
  -w0 * (-muPsi / sigmaPsi * sqrt(2 / pi) * exp(-0.5 * muPsi^2 / sigmaPsi2) + (1 - 2 * pnorm(-muPsi / sigmaPsi)) + 2 * muPsi / sigmaPsi * dnorm(-muPsi / sigmaPsi))
}

DS2_muPsi <- function(muPsi, sigmaPsi2, sigmaPsi, w0, rqt, sqt, SigmaThetaq, muThetaq, J) {
  -0.25 * J * (J + 1) / w0 * DS1_sigmaPsi2(muPsi, sigmaPsi2, sigmaPsi, w0) - 0.5 * rqt / sqt * sum((diag(SigmaThetaq) + muThetaq^2) * DQj_muPsi(muPsi, sigmaPsi2, sigmaPsi, 1:J))
}

DS2_sigmaPsi2 <- function(muPsi, sigmaPsi2, sigmaPsi, w0, rqt, sqt, SigmaThetaq, muThetaq, J) {
  -0.25 * J * (J + 1) / w0 * DS1_sigmaPsi2(muPsi, sigmaPsi2, sigmaPsi, w0) - 0.5 * rqt / sqt * sum((diag(SigmaThetaq) + muThetaq^2) * DQj_sigmaPsi2(muPsi, sigmaPsi2, sigmaPsi, 1:J))
}

# LB <- function(y, W, varphi, WtW, varphitvarphi, muPsi, sigmaPsi2, sigmaPsi, w0, J, SigmaThetaq, muThetaq, SigmaBetaq, muBetaq, SigmaBeta0_inv, muBeta0, rqt, sqt, rqt_half, sqt_half, rqs, sqs, rqs_half, sqs_half, r0t_half, s0t_half, r0s_half, s0s_half) {
#   p <- dim(W)[2]
#   t1 <- SigmaBeta0_inv %*% SigmaBetaq
#   -0.5 * (sum(WtW * SigmaBetaq) + sum(varphitvarphi * SigmaThetaq)) + sum(log(pnorm(W %*% muBetaq + varphi %*% muThetaq)^y * (1 - pnorm(W %*% muBetaq + varphi %*% muThetaq))^(1 - y))) - 0.5 * J * (log(2 * pi) - (digamma(rqs_half) - log(sqs_half)) - (digamma(rqt_half) - log(sqt_half))) + S2(muPsi, sigmaPsi2, sigmaPsi, w0, rqs, sqs, rqt, sqt, SigmaThetaq, muThetaq, J) + 0.5 * J * (1 + log(2 * pi)) + 0.5 * determinant(SigmaThetaq)$modulus[1] + r0t_half * log(r0t_half) - lgamma(r0t_half) + (r0t_half - rqt_half) * digamma(rqt_half) - (r0t_half + 1) * log(sqt_half) - r0t_half * rqt / sqt + rqt_half + log(sqt_half) + lgamma(rqt_half) + r0s_half * log(s0s_half) - lgamma(r0s_half) + (r0s_half - rqs_half) * digamma(rqs_half) - (r0s_half + 1) * log(sqs_half) - sqs_half * rqs / sqs + rqs_half + log(sqs_half) + lgamma(rqs_half) + 0.5 * (p + 1) * (1 + digamma(rqs_half) - log(sqs_half)) + 0.5 * determinant(t1)$modulus[1] - 0.5 * rqs / sqs * (sum(diag(t1)) + sum((muBetaq - muBeta0) * (SigmaBeta0_inv %*% (muBetaq - muBeta0)))) + log(w0 / 2) + S1(muPsi, sigmaPsi2, sigmaPsi, w0) + 0.5 * (log(2 * pi * sigmaPsi2) - 1)
# }


LB <- function(y, muYStar, W, varphi, WtW, varphitvarphi, muPsi, sigmaPsi2, sigmaPsi, w0, J, SigmaThetaq, muThetaq, SigmaBetaq, muBetaq, SigmaBeta0_inv, muBeta0, rqt, sqt, rqt_half, sqt_half, r0t_half, s0t_half) {
  p <- dim(W)[2]
  t1 <- SigmaBeta0_inv %*% SigmaBetaq
  -0.5 * (sum(WtW * SigmaBetaq) + sum(varphitvarphi * SigmaThetaq)) + sum(muYStar * dnorm(muYStar) / (2 * pnorm(muYStar)^y * (pnorm(muYStar) - 1)^(1 - y))) - sum(muYStar * dnorm(muYStar) / (2 * pnorm(muYStar)^y * (1 - pnorm(muYStar))^(1 - y))) + 0.5 * (1 + log(2 * pi)) + sum(log(pnorm(muYStar)^y * (1 - pnorm(muYStar))^(1 - y)))  - 0.5 * J * (log(2 * pi) - (digamma(rqt_half) - log(sqt_half))) + S2(muPsi, sigmaPsi2, sigmaPsi, w0, rqt, sqt, SigmaThetaq, muThetaq, J) + 0.5 * J * (1 + log(2 * pi)) + 0.5 * determinant(SigmaThetaq)$modulus[1] + r0t_half * log(r0t_half) - lgamma(r0t_half) + (r0t_half - rqt_half) * digamma(rqt_half) - (r0t_half + 1) * log(sqt_half) - r0t_half * rqt / sqt + rqt_half + log(sqt_half) + lgamma(rqt_half) + 0.5 * determinant(t1)$modulus[1] - 0.5 * (sum(diag(t1)) + sum((muBetaq - muBeta0) * (SigmaBeta0_inv %*% (muBetaq - muBeta0)))) + log(w0 / 2) + S1(muPsi, sigmaPsi2, sigmaPsi, w0) + 0.5 * (log(2 * pi * sigmaPsi2) - 1)
  # -0.5 * (sum(WtW * SigmaBetaq) + sum(varphitvarphi * SigmaThetaq)) + sum(log(pnorm(W %*% muBetaq + varphi %*% muThetaq)^y * (1 - pnorm(W %*% muBetaq + varphi %*% muThetaq))^(1 - y))) - 0.5 * J * (log(2 * pi) - (digamma(rqs_half) - log(sqs_half)) - (digamma(rqt_half) - log(sqt_half))) + S2(muPsi, sigmaPsi2, sigmaPsi, w0, rqs, sqs, rqt, sqt, SigmaThetaq, muThetaq, J) + 0.5 * J * (1 + log(2 * pi)) + 0.5 * determinant(SigmaThetaq)$modulus[1] + r0t_half * log(r0t_half) - lgamma(r0t_half) + (r0t_half - rqt_half) * digamma(rqt_half) - (r0t_half + 1) * log(sqt_half) - r0t_half * rqt / sqt + rqt_half + log(sqt_half) + lgamma(rqt_half) + r0s_half * log(s0s_half) - lgamma(r0s_half) + (r0s_half - rqs_half) * digamma(rqs_half) - (r0s_half + 1) * log(sqs_half) - sqs_half * rqs / sqs + rqs_half + log(sqs_half) + lgamma(rqs_half) + 0.5 * (p + 1) * (1 + digamma(rqs_half) - log(sqs_half)) + 0.5 * determinant(t1)$modulus[1] - 0.5 * rqs / sqs * (sum(diag(t1)) + sum((muBetaq - muBeta0) * (SigmaBeta0_inv %*% (muBetaq - muBeta0)))) + log(w0 / 2) + S1(muPsi, sigmaPsi2, sigmaPsi, w0) + 0.5 * (log(2 * pi * sigmaPsi2) - 1)
}

# LB <- function(y, W, varphi, WtW, varphitvarphi, muPsi, sigmaPsi2, sigmaPsi, w0, J, SigmaThetaq, muThetaq, SigmaBetaq, muBetaq, SigmaBeta0_inv, muBeta0, muYStar, muYStarq, rqt, sqt, rqt_half, sqt_half, rqs, sqs, rqs_half, sqs_half, r0t_half, s0t_half, r0s_half, s0s_half) {
# p <- dim(W)[2]
# t1 <- SigmaBeta0_inv %*% SigmaBetaq
#   -0.5 * (sum(muYStar^2) + 1 + sum(muYStar * (dnorm(muYStar) / ((pnorm(muYStar)^y) * (pnorm(muYStar) - 1)^(1 - y)))) - 2 * sum(muYStarq * (W %*% muBetaq + varphi %*% muThetaq)) + sum((W %*% muBetaq + varphi %*% muThetaq)^2) + sum(WtW * SigmaBetaq) + sum(varphitvarphi * SigmaThetaq)) + 0.5 * (1 + log(2 * pi)) + sum(log(pnorm(muYStar)^y * (1 - pnorm(muYStar)^(1 - y)))) - sum(muYStar * dnorm(muYStar) / (2 * pnorm(muYStar)^y * (1 - pnorm(muYStar))^(1 - y))) - 0.5 * J * (log(2 * pi) - (digamma(rqs_half) - log(sqs_half)) - (digamma(rqt_half) - log(sqt_half))) + S2(muPsi, sigmaPsi2, sigmaPsi, w0, rqs, sqs, rqt, sqt, SigmaThetaq, muThetaq, J) + 0.5 * J * (1 + log(2 * pi)) + 0.5 * determinant(SigmaThetaq)$modulus[1] + r0t_half * log(r0t_half) - lgamma(r0t_half) + (r0t_half - rqt_half) * digamma(rqt_half) - (r0t_half + 1) * log(sqt_half) - r0t_half * rqt / sqt + rqt_half + log(sqt_half) + lgamma(rqt_half) + r0s_half * log(s0s_half) - lgamma(r0s_half) + (r0s_half - rqs_half) * digamma(rqs_half) - (r0s_half + 1) * log(sqs_half) - sqs_half * rqs / sqs + rqs_half + log(sqs_half) + lgamma(rqs_half) + 0.5 * (p + 1) * (1 + digamma(rqs_half) - log(sqs_half)) + 0.5 * determinant(t1)$modulus[1] - 0.5 * rqs / sqs * (sum(diag(t1)) + sum((muBetaq - muBeta0) * (SigmaBeta0_inv %*% (muBetaq - muBeta0)))) + log(w0 / 2) + S1(muPsi, sigmaPsi2, sigmaPsi, w0) + 0.5 * (log(2 * pi * sigmaPsi2) - 1)
# }


# # Lower bound before NCVMP
# LBR <- function() {
#   t1 <- SigmaBeta0_inv %*% SigmaBetaq
#   -0.5 * (sum(WtW * SigmaBetaq) + sum(varphitvarphi * SigmaThetaq)) + sum(log(pnorm(W %*% muBetaq + varphi %*% muThetaq)^y * (1 - pnorm(W %*% muBetaq + varphi %*% muThetaq))^(1 - y))) - 0.5 * J * (log(2 * pi) - (digamma(rqs_half) - log(sqs_half)) - (digamma(rqt_half) - log(sqt_half))) + 0.5 * J * (1 + log(2 * pi) + determinant(SigmaThetaq)$modulus[1]) + r0t_half * log(s0t_half) - lgamma(r0t_half) + (r0t_half - rqt_half) * digamma(rqt_half) - (r0t_half + 1) * log(sqt_half) - r0s_half * log(s0s_half) - lgamma(r0s_half) + (r0s_half - rqs_half) * digamma(rqs_half) - (r0s_half + 1) * log(sqs_half) - s0s_half * rqs / sqs + rqs_half + log(sqs_half) + lgamma(rqs_half) + 0.5 * p * (1 + digamma(rqs_half) - log(sqs_half)) + 0.5 * determinant(t1)$modulus[1] - 0.1
# }


VB <- function(x, y, W, muBeta0, SigmaBeta0, w0, r0t, s0t, J, tol = 1.0e-06) {
  n <- length(y)
  if (!is.matrix(W)) W <- as.matrix(W)
  p <- dim(W)[2]
  dif <- tol + 1
  varphi <- sqrt(2) * cos(outer(x, pi * (1:J)))
  varphitvarphi <- crossprod(varphi)
  WtW <- crossprod(W)
  SigmaBeta0_inv <- solve(SigmaBeta0)
  SigmaBeta0_inv_muBeta0 <- solve(SigmaBeta0, muBeta0)
  # rqs <- r0s + J + p
  # rqs_half <- 0.5 * rqs
  # sqs <- s0s + 200
  # sqs_half <- 0.5 * sqs
  rqt <- r0t + J
  rqt_half <- 0.5 * rqt
  sqt <- s0t + 200
  sqt_half <- 0.5 * sqt
  s0t_half <- 0.5 * s0t
  # s0s_half <- 0.5 * s0s
  r0t_half <- 0.5 * r0t
  # r0s_half <- 0.5 * r0s
  muBetaq <- muBeta0
  muPsi <- 0.01
  sigmaPsi2 <- muPsi^2 / 100
  sigmaPsi <- sqrt(sigmaPsi2)
  muThetaq <- rep(0, J)
  muYStar <- rep(0, n)
  muYStarq <- rep(0.1, n)
  SigmaBetaq <- SigmaBeta0
  SigmaThetaq <- diag(1, J)
  count <- 0
  
  lbold <- LB(y, muYStar, W, varphi, WtW, varphitvarphi, muPsi, sigmaPsi2, sigmaPsi, w0, J, SigmaThetaq, muThetaq, SigmaBetaq, muBetaq, SigmaBeta0_inv, muBeta0, rqt, sqt, rqt_half, sqt_half, r0t_half, s0t_half)
  # lbold <- LB(y, W, varphi, WtW, varphitvarphi, muPsi, sigmaPsi2, sigmaPsi, w0, J, SigmaThetaq, muThetaq, SigmaBetaq, muBetaq, SigmaBeta0_inv, muBeta0, rqt, sqt, rqt_half, sqt_half, rqs, sqs, rqs_half, sqs_half, r0t_half, s0t_half, r0s_half, s0s_half)
  cat('first lb: ', lbold, '\n')
  lbnew <- 0
  lbrecord <- c(lbold)
  while (dif > tol | dif < 0) {
    count <- count + 1

    # Update theta
    # cat('rqs / sqs: ', rqs / sqs, '\n')
    # cat('rqt / sqt: ', rqt / sqt, '\n')
    # cat('Qj(muPsi, sigmaPsi2, sigmaPsi, 1:J): ', Qj(muPsi, sigmaPsi2, sigmaPsi, 1:J), '\n')
    # SigmaThetaq <- solve(varphitvarphi + rqs / sqs * rqt / sqt * diag(Qj(muPsi, sigmaPsi2, sigmaPsi, 1:J)))
    SigmaThetaq <- solve(varphitvarphi + rqt / sqt * diag(Qj(muPsi, sigmaPsi2, sigmaPsi, 1:J)))
    muThetaq <- SigmaThetaq %*% crossprod(varphi, muYStarq - W %*% muBetaq)

    # Update tau
    sqt <- s0t + sum((diag(SigmaThetaq) + muThetaq^2) * Qj(muPsi, sigmaPsi2, sigmaPsi, 1:J))
    sqt_half <- 0.5 * sqt

    # Update sigma
    # cat('sum((diag(SigmaThetaq) + muThetaq^2) * Qj(muPsi, sigmaPsi2, sigmaPsi, 1:J)):', sum((diag(SigmaThetaq) + muThetaq^2) * Qj(muPsi, sigmaPsi2, sigmaPsi, 1:J)), '\n')
    # cat('sum(SigmaBeta0_inv * SigmaBetaq): ', sum(SigmaBeta0_inv * SigmaBetaq), '\n')
    # sqs <- s0s + rqt / sqt * sum((diag(SigmaThetaq) + muThetaq^2) * Qj(muPsi, sigmaPsi2, sigmaPsi, 1:J)) + sum(SigmaBeta0_inv * SigmaBetaq) + sum((muBetaq - muBeta0) * (SigmaBeta0_inv %*% (muBetaq - muBeta0)))

    # Update beta
    # cat('rqs/sqs: ', rqs/sqs, '\n')
    # cat('SigmaBeta0_inv: ', SigmaBeta0_inv, '\n')
    # cat('WtW: ', WtW, '\n')
    # cat('rqs/sqs * SigmaBeta0_inv + WtW: ', rqs / sqs * SigmaBeta0_inv + WtW, '\n')
    SigmaBetaq <- solve(as.matrix(SigmaBeta0_inv + WtW))
    muBetaq <- SigmaBetaq %*% (SigmaBeta0_inv_muBeta0 + crossprod(muYStar - varphi %*% muThetaq))

    # Update muYStar
    muYStar <- W %*% muBetaq + varphi %*% muThetaq
    # for (i in 1:n) {
    #   if (y[i] == 1 & muYStar[i] < -8) {
    #     muYStarq[i] <- muYStar[i] + dnormOverPnorm_Laurent(muYStar[i])
    #   } else if (y[i] == 1 & muYStar[i] >= -8) {
    #     muYStarq[i] <- muYStar[i] + dnormOverPnorm(muYStar[i])
    #   } else {
    #     muYStarq[i] <- muYStar[i] + dnormOverPnormMinusOne(muYStar[i])
    #   }
    # }
    muYStarq <- muYStar + dnorm(muYStar)/((pnorm(muYStar)^y)*(pnorm(muYStar)-1)^(1-y))

    ## Calculate lower bound before NCVMP
    lbfull <- LB(y, muYStar, W, varphi, WtW, varphitvarphi, muPsi, sigmaPsi2, sigmaPsi, w0, J, SigmaThetaq, muThetaq, SigmaBetaq, muBetaq, SigmaBeta0_inv, muBeta0, rqt, sqt, rqt_half, sqt_half, r0t_half, s0t_half)
    # lbfull <- LB(y, W, varphi, WtW, varphitvarphi, muPsi, sigmaPsi2, sigmaPsi, w0, J, SigmaThetaq, muThetaq, SigmaBetaq, muBetaq, SigmaBeta0_inv, muBeta0, rqt, sqt, rqt_half, sqt_half, rqs, sqs, rqs_half, sqs_half, r0t_half, s0t_half, r0s_half, s0s_half)

    # Update psi
    a <- 1
    sigmaPsi2_old <- sigmaPsi2
    muPsi_old <- muPsi
    temp <- DS1_muPsi(muPsi, sigmaPsi2, sigmaPsi, w0) + DS2_muPsi(muPsi, sigmaPsi2, sigmaPsi, w0, rqt, sqt, SigmaThetaq, muThetaq, J)
    # sigmaPsi2 <- -0.5 / (DS1_sigmaPsi2(muPsi, sigmaPsi2, sigmaPsi, w0) + DS2_sigmaPsi2(muPsi, sigmaPsi2, sigmaPsi, w0, rqs, sqs, rqt, sqt, SigmaThetaq, muThetaq, J))
    
    sigmaPsi2_new <- -0.5 / (DS1_sigmaPsi2(muPsi, sigmaPsi2, sigmaPsi, w0) + DS2_sigmaPsi2(muPsi, sigmaPsi2, sigmaPsi, w0, rqt, sqt, SigmaThetaq, muThetaq, J))
    # temp2 <- sigmaPsi2_new

    while (sigmaPsi2_new < 0) {
      cat('sigmaPsi2 is negative! \n')
      a <- 2/3 * a
      sigmaPsi2_new <- 1 / (1 / sigmaPsi2_old + a * (1 / sigmaPsi2_new - 1 / sigmaPsi2_old))
    }
    sigmaPsi2 <- sigmaPsi2_new
    sigmaPsi <- sqrt(sigmaPsi2)
    muPsi <- muPsi + a * sigmaPsi2 * temp
    # cat('sigmaPsi2: ', sigmaPsi2, '\n')
    # cat('muPsi: ', muPsi, '\n')
    lbnew <- LB(y, muYStar, W, varphi, WtW, varphitvarphi, muPsi, sigmaPsi2, sigmaPsi, w0, J, SigmaThetaq, muThetaq, SigmaBetaq, muBetaq, SigmaBeta0_inv, muBeta0, rqt, sqt, rqt_half, sqt_half, r0t_half, s0t_half)
    # lbnew <- LB(y, W, varphi, WtW, varphitvarphi, muPsi, sigmaPsi2, sigmaPsi, w0, J, SigmaThetaq, muThetaq, SigmaBetaq, muBetaq, SigmaBeta0_inv, muBeta0, rqt, sqt, rqt_half, sqt_half, rqs, sqs, rqs_half, sqs_half, r0t_half, s0t_half, r0s_half, s0s_half)
    dif <- lbnew - lbold
    cat('lbfull: ', lbfull, '\n')

    if (dif < 0) {
      a <- 1
      dif_try <- dif
      while (dif_try < 0) {
        a <- 2/3 * a
        sigmaPsi2_try <- 1 / (1 / sigmaPsi2_old + a * (1 / sigmaPsi2 - 1 / sigmaPsi2_old))
        sigmaPsi_try <- sqrt(sigmaPsi2_try)
        muPsi_try <- sigmaPsi2_try * (muPsi_old / sigmaPsi2_old + a * (muPsi / sigmaPsi2 - muPsi_old / sigmaPsi2_old))
        lbnew <- LB(y, muYStar, W, varphi, WtW, varphitvarphi, muPsi_try, sigmaPsi2_try, sigmaPsi_try, w0, J, SigmaThetaq, muThetaq, SigmaBetaq, muBetaq, SigmaBeta0_inv, muBeta0, rqt, sqt, rqt_half, sqt_half, r0t_half, s0t_half)
        # lbnew <- LB(y, W, varphi, WtW, varphitvarphi, muPsi_try, sigmaPsi2_try, sigmaPsi_try, w0, J, SigmaThetaq, muThetaq, SigmaBetaq, muBetaq, SigmaBeta0_inv, muBeta0, rqt, sqt, rqt_half, sqt_half, rqs, sqs, rqs_half, sqs_half, r0t_half, s0t_half, r0s_half, s0s_half)
        dif_try <- lbnew - lbfull
      }

      sigmaPsi2 <- sigmaPsi2_try
      sigmaPsi <- sigmaPsi_try
      muPsi <- muPsi_try
    #   # Update psi
    #   sigmaPsi2_old <- sigmaPsi2
    #   muPsi_old <- muPsi
    #   sigmaPsi2_new <- -0.5 / (DS1_sigmaPsi2(muPsi, sigmaPsi2, sigmaPsi, w0) + DS2_sigmaPsi2(muPsi, sigmaPsi2, sigmaPsi, w0, rqs, sqs, rqt, sqt, SigmaThetaq, muThetaq, J))
    #   temp2 <- sigmaPsi2_new

    #   temp <- DS1_muPsi(muPsi, sigmaPsi2, sigmaPsi, w0) + DS2_muPsi(muPsi, sigmaPsi2, sigmaPsi, w0, rqs, sqs, rqt, sqt, SigmaThetaq, muThetaq, J)
    #   while (sigmaPsi2_new < 0) {
    #     a <- 2/3 * a
    #     sigmaPsi2_new <- 1 / (1 / sigmaPsi2_old + a * (1 / sigmaPsi2_new - 1 / sigmaPsi2_old))
    #   }
    #   sigmaPsi2 <- sigmaPsi2_new
    #   sigmaPsi <- sqrt(sigmaPsi2)
    #   muPsi <- muPsi + a * sigmaPsi2 * temp

    #   # Update theta
    #   SigmaThetaq <- solve(varphitvarphi + rqs / sqs * rqt / sqt * diag(Qj(muPsi, sigmaPsi2, sigmaPsi, 1:J)))
    #   muThetaq <- SigmaThetaq %*% crossprod(varphi, muYStarq - W %*% muBetaq)

    #   # Update tau
    #   sqt <- s0t + rqs / sqs * sum((diag(SigmaThetaq) + muThetaq^2) * Qj(muPsi, sigmaPsi2, sigmaPsi, 1:J))
    #   sqt_half <- 0.5 * sqt

    #   # Update sigma
    #   sqs <- s0s + rqt / sqt * sum((diag(SigmaThetaq) + muThetaq^2) * Qj(muPsi, sigmaPsi2, sigmaPsi, 1:J)) + sum(SigmaBeta0_inv * SigmaBetaq) + sum((muBetaq - muBeta0) * (SigmaBeta0_inv %*% (muBetaq - muBeta0)))

    #   # Update beta
    #   SigmaBetaq <- solve(rqs / sqs * SigmaBeta0_inv + WtW)
    #   muBetaq <- SigmaBetaq %*% (rqs / sqs * SigmaBeta0_inv_muBeta0 + crossprod(muYStar - varphi %*% muThetaq))

    #   # Update muYStar
    #   muYStar <- W %*% muBetaq + varphi %*% muThetaq
    #   for (i in 1:n) {
    #     if (y[i] > 0 & muYStar[i] < -8) {
    #       muYStarq[i] <- muYStar[i] + dnormOverPnorm_Laurent(muYStar[i])
    #     } else if (y[i] > 0 & muYStar[i] >= -8) {
    #       muYStarq[i] <- muYStar[i] + dnormOverPnorm(muYStar[i])
    #     } else {
    #       muYStarq[i] <- muYStar[i] + dnormOverPnormMinusOne(muYStar[i])
    #     }
    #   }
    #   lbnew <- LB(y, W, varphi, WtW, varphitvarphi, muPsi, sigmaPsi2, sigmaPsi, w0, J, SigmaThetaq, muThetaq, SigmaBetaq, muBetaq, SigmaBeta0_inv, muBeta0, rqt, sqt, rqt_half, sqt_half, rqs, sqs, rqs_half, sqs_half, r0t_half, s0t_half, r0s_half, s0s_half)
    #   diff <- lbnew - lbold
    # }

    # dif <- (lbnew - lbold) / abs(lbnew)
    # lbold <- lbnew
    # lbrecord <- c(lbrecord, lbnew)
    }
    dif <- (lbnew - lbold)/abs(lbnew)
    # dif <- abs(lbnew - lbold)
    lbold <- lbnew
    lbrecord <- c(lbrecord, lbnew)
    cat('count: ', count, ', lbnew: ', lbnew, ', dif: ', dif, '\n')
  }
  list(muThetaq = muThetaq, SigmaThetaq = SigmaThetaq, muBetaq = muBetaq, SigmaBetaq = SigmaBetaq, rqt = rqt, sqt = sqt, muYStarq = muYStarq, sigmaPsi2 = sigmaPsi2, muPsi = muPsi, lbrecord = lbrecord, count = count)
}