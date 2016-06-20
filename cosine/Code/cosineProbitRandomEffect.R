# Auxiliary functions
Qj = function(muPsi, sigmaPsi2, j) {
  sigmaPsi = sqrt(sigmaPsi2)
  exp(0.5 * sigmaPsi2 * j^2 + muPsi * j) * (1 - pnorm(-muPsi/sigmaPsi - sigmaPsi * j)) + exp(9.5 * sigmaPsi2 * j^2 - muPsi * j) * (1 - pnorm(muPsi/sigmaPsi - sigmaPsi * j))
}

S1 = function(muPsi, sigmaPsi2, w0) {
  sigmaPsi = sqrt(sigmaPsi2)
  -w0 * (sigmaPsi * sqrt(2 / pi) * exp(-muPsi^2 / (2 * sigmaPsi2)) + muPsi * (1 - 2 * pnorm(-muPsi / sigmaPsi)))
}

S2 = function(muPsi, sigmaPsi2, muThetaq, SigmaThetaq, rqs, sqs, rqt, sqt, w0, J) {
  -0.5 * rqs/sqs * rqt/sqt * sum((diag(SigmaThetaq) + muThetaq^2) * Qj(muPsi, sigmaPsi2, 1:J)) - 0.25 * J * (J + 1) / w0 * S1(muPsi, sigmaPsi2, w0)
}

DS1_sigmaPsi2 = function(muPsi, sigmaPsi2, w0) {
  sigmaPsi = sqrt(sigmaPsi2)
  -w0 * ((1/(2 * sigmaPsi) + 0.5 * sigmaPsi * muPsi^2 / sigmaPsi2^2) * sqrt(2 / pi) * exp(-0.5 * muPsi^2 / sigmaPsi2) - muPsi^2/sigmaPsi^3 * dnorm(-muPsi/sigmaPsi))
}

DS1_muPsi = function(muPsi, sigmaPsi2, w0) {
  sigmaPsi = sqrt(sigmaPsi2)
  -w0 * (-muPsi/sigmaPsi * sqrt(2/pi) * exp(-0.5 * muPsi^2 / sigmaPsi2) + (1 - 2 * pnorm(-muPsi/sigmaPsi)) + 2 * muPsi / sigmaPsi * dnorm(-muPsi/sigmaPsi))
}

DQj_muPsi = function(muPsi, sigmaPsi2, j) {
  sigmaPsi = sqrt(sigmaPsi2)
  term1 = 0.5 * sigmaPsi2 * j^2
  term2 = muPsi * j
  term3 = muPsi / sigmaPsi
  term4 = sigmaPsi * j
  term5 = 1/sigmaPsi
  exp(term1 + term2) * dnorm(-term3 - term4) * term5 + j * exp(term1 + term2) * (1 - pnorm(-term3 - term4)) - exp(term1 - term2) * dnorm(term3 - term4) * term5 - j * exp(term1 - term2) * (1 - pnorm(term3 - term4))
}

DQj_sigmaPsi2 = function(muPsi, sigmaPsi2, j) {
  sigmaPsi = sqrt(sigmaPsi2)
  term1 = 0.5 * sigmaPsi2 * j^2
  term2 = muPsi * j
  term3 = muPsi / sigmaPsi
  term4 = sigmaPsi * j
  term5 = 0.5 * muPsi / sigmaPsi^3
  term6 = 0.5 * j / sigmaPsi2
  term7 = 0.5 * j^2
  -exp(term1 + term2) * dnorm(-term3 - term4) * (term5 - term6) + term7 * exp(term1 + term2) * (1 - pnorm(-term3 - term4)) + exp(term1 - term2) * dnorm(term3 - term4) * (term5 + term6) + term7 * exp(term1 - term2) * (1 - pnorm(term3 - term4))
}

DS2_sigmaPsi2 = function(muPsi, sigmaPsi2, w0, rqs_over_sqs, rqt_over_sqt, SigmaThetaq, muThetaq) {
  J = length(muThetaq)
  -0.25 * J * (J + 1) / w0 * DS1_sigmaPsi2(muPsi, sigmaPsi2, w0) - 0.5 * rqs_over_sqs * rqt_over_sqt * sum((diag(SigmaThetaq) + muThetaq^2) * DQj_sigmaPsi2(muPsi, sigmaPsi2, 1:J))
}

DS2_muPsi = function(muPsi, sigmaPsi2, w0, rqs_over_sqs, rqt_over_sqt, SigmaThetaq, muThetaq) {
  J = length(muThetaq)
  -0.25 * J * (J + 1) / w0 * DS1_muPsi(muPsi, sigmaPsi2, w0) - 0.5 * rqs_over_sqs * rqt_over_sqt * sum((diag(SigmaThetaq) + muThetaq^2) * DQj_muPsi(muPsi, sigmaPsi2, 1:J))
}

LB = function(y, W, Z, varphi, SigmaBetaq, Sigmauq, SigmaThetaq, muYStar, Sigmau0, rqs, sqs, rqt, sqt, muPsi, sigmaPsi2, r0t, s0t, r0s, s0s, , SigmaBeta0, muBetaq, muBeta0, w0, muThetaq) {
  J = length(muThetaq)
  p = length(muBetaq)
  t1 = solve(SigmaBeta0, SigmaBetaq)
  -0.5 * (sum(crossprod(W) * SigmaBetaq) + sum(crossprod(Z) * Sigmauq) + sum(crossprod(varphi) * SigmaThetaq)) + sum(log((pnorm(muYStar))^y * (1 - pnorm(muYStar))^(1-y))) + 0.5 * determinant(solve(Sigmau0, Sigmauq))$modulus[1] - 0.5 * J * (log(2 * pi) - (digamma(rqs/2) - log(sqs/2)) - (digamma(rqt/2) - log(sqt/2))) + S2(muPsi, sigmaPsi2, muThetaq, SigmaThetaq, rqs, sqs, rqt, sqt, w0, J) + r0t/2 * log(s0t/2) - lgamma(r0t/2) + (r0t/2 + 1) * (digamma(rqt/2) - log(sqt/2)) - s0t/2 * rqt/sqt + rqt/2 + log(sqt/2) + lgamma(rqt/2) - (1 + rqt/2) * digamma(rqt/2) + r0s/2 * log(s0s/2) - lgamma(r0s/2) + (r0s/2 + 1) * (digamma(rqs/2) - log(sqs/2)) - s0s/2 * rqs/sqs + rqs/2 + log(sqs/2) + lgamma(rqs/2) - (1 + rqs/2) * digamma(rqs/2) + 0.5 * (p + 1) + 0.5 * (digamma(rqs/2) - log(sqs/2)) + 0.5 * determinant(t1)$modulus[1] - 0.5 * rqs/sqs * (sum(diag(t1)) + sum((muBetaq - muBeta0) * solve(SigmaBeta0, muBetaq - muBeta0))) + log(w0/2) + S1(muPsi, sigmaPsi2, w0) + 0.5 * log(2 * pi * sigmaPsi2) - 0.5
}

VB = function(y, x, W, muBeta0, SigmaBeta0, w0, r0s, s0s, r0t, s0t, muPsi0, muu0, Sigmau0, J, tol = 1.0e-06) {
  n = length(y)
  p = dim(W)[2]
  dif = tol + 1
  varphi = sqrt(2)*cos(outer(x,pi*(1:J)))
  varphitvarphi = crossprod(varphi)
  WtW = crossprod(W)
  SigmaBeta0_inv = solve(SigmaBeta0)
  SigmaBeta0_inv_muBeta0 = SigmaBeta0_inv %*% muBeta0
  r0t_over_s0t = r0t / s0t
  r0s_over_s0s = r0s / s0s
  r0t_half = 0.5 * r0t
  s0t_half = 0.5 * s0t
  r0s_half = 0.5 * r0s
  s0s_half = 0.5 * s0s
  # Initialize variational parameters
  rqs = r0s + J + p + 1
  sqs = r0s
  rqt = r0t + J
  sqt = s0t
  rqt_over_sqt = rqt / sqt
  rqs_over_sqs = rqs / sqs
  muuq = muu0
  Sigmauq = Sigmau0
  muBetaq = muBeta0
  muPsi = muPsi0
  sigmaPsi2 = muPsi^2 / 10
  sigmaPsi = sqrt(sigmaPsi2)
  muThetaq = rnorm(J, 0, s0t_half / (r0t_half - 1) * s0s_half / (r0s_half - 1) * exp(-(1:J) * abs(muPsi0)))
  muYStar = W %*% muBetaq + Z %*% muuq + varphi %*% muThetaq
  muYStarq = muYStar + dnorm(muYStar) / (pnorm(muYStar)^y * (pnorm(muYStar) - 1)^(1-y))
  SigmaThetaq = solve(varphitvarphi + rqs_over_sqs * rqt_over_sqt * Qj(muPsi, sigmaPsi2, 1:J))

  count = 0
  lbold = LB(y, W, Z, varphi, SigmaBetaq, Sigmauq, SigmaThetaq, muYStar, Sigmau0, rqs, sqs, rqt, sqt, muPsi, sigmaPsi2, r0t, s0t, r0s, s0s, , SigmaBeta0, muBetaq, muBeta0, w0, muThetaq)
  lbrecord = c(lbold)
  while (dif > tol) {
    count = count + 1
    a = 1
    # Update psi
    sigmaPsi2_old = sigmaPsi2
    muPsi_old = muPsi
    sigmaPsi2_new = -0.5 / (DS1_sigmaPsi2(muPsi, sigmaPsi2, w0) + DS2_sigmaPsi2(muPsi, sigmaPsi2, w0, rqs_over_sqs, rqt_over_sqt, SigmaThetaq, muThetaq))
    temp2 = sigmaPsi2_new
    muPsi_old = muPsi

    temp = DS1_muPsi(muPsi, sigmaPsi2, w0) + DS2_muPsi(muPsi, sigmaPsi2, w0, rqs_over_sqs, rqt_over_sqt, SigmaThetaq, muThetaq)
    while (sigmPsi2_new < 0) {
      a = 2/3 * a
      sigmaPsi2_new = 1 / (1/sigmaPsi2_old + a * (1/sigmaPsi2_new - 1/sigmaPsi2_old))
    }
    sigmaPsi2 = sigmaPsi2_new
    sigmaPsi = sqrt(sigmaPsi2)
    muPsi = muPsi + a * sigmaPsi2 * temp

    # Update theta
    SigmaThetaq = solve(varphitvarphi + (rqs_over_sqs * rqt_over_sqt * Qj(muPsi, sigmaPsi2, 1:J)))
    muThetaq = SigmaThetaq %*% crossprod(varphi, muYStarq - W %*% muThetaq - Z %*% muuq)

    # Update tau^2
    sqt = s0t + (rqt_over_sqt * sum((diag(SigmaThetaq) + muThetaq^2) * Qj(muPsi, sigmaPsi2, 1:J)))
    sqt_half = 0.5 * sqt
    rqt_over_sqt = rqt / sqt

    # Update sigma^2
    sqs = s0s + rqt_over_sqt * sum((diag(SigmaThetaq) + muThetaq^2) * Qj(muPsi, sigmaPsi2, 1:J)) + sum(diag(solve(SigmaBeta0, SigmaBetaq))) + sum(SigmaBeta0_inv * SigmaBetaq) + sum(SigmaBeta0_inv * tcrossprod(muBetaq - muBeta0))
    sqs_half = 0.5 * sqs
    rqs_over_sqs = rqs / sqs

    # Update beta
    SigmaBetaq = solve(rqs_over_sqs * SigmaBeta0_inv + WtW)
    muBetaq = SigmaBetaq %*% (SigmaBeta0_inv_muBeta0 + crossprod(W, muYStarq - Z %*% muuq - varphi %*% muThetaq))

    # Update muYStar
    muYStar = W %*% muBetaq + Z %*% muuq + varphi %*% muThetaq
    muYStarq = muYStar + dnorm(muYStar) / (pnorm(muYStar)^y * (pnorm(muYStar) - 1)^(1-y))

    lbnew = LB(y, W, Z, varphi, SigmaBetaq, Sigmauq, SigmaThetaq, muYStar, Sigmau0, rqs, sqs, rqt, sqt, muPsi, sigmaPsi2, r0t, s0t, r0s, s0s, , SigmaBeta0, muBetaq, muBeta0, w0, muThetaq)
    diff = lbnew - lbold
    if (diff < 0) {
      a = 1
      sigmPsi2_new = temp2
      while (sigmaPsi_new < 0) {
        a = 2/3 * a
        sigmaPsi_new = 1/(1/sigmaPsi2 + a * (1/temp2 - 1/sigmaPsi2))
      }
      sigmaPsi2 = sigmaPsi2_new
      sigmaPsi = sqrt(sigmaPsi2)
      muPsi = muPsi_old + a * sigmaPsi2 * temp

      # Update theta
      SigmaThetaq = solve(varphitvarphi + (rqs_over_sqs * rqt_over_sqt * Qj(muPsi, sigmaPsi2, 1:J)))
      muThetaq = SigmaThetaq %*% crossprod(varphi, muYStarq - W %*% muThetaq - Z %*% muuq)

      # Update tau^2
      sqt = s0t + (rqt_over_sqt * sum((diag(SigmaThetaq) + muThetaq^2) * Qj(muPsi, sigmaPsi2, 1:J)))
      sqt_half = 0.5 * sqt
      rqt_over_sqt = rqt / sqt

      # Update sigma^2
      sqs = s0s + rqt_over_sqt * sum((diag(SigmaThetaq) + muThetaq^2) * Qj(muPsi, sigmaPsi2, 1:J)) + sum(diag(solve(SigmaBeta0, SigmaBetaq))) + sum(SigmaBeta0_inv * SigmaBetaq) + sum(SigmaBeta0_inv * tcrossprod(muBetaq - muBeta0))
      sqs_half = 0.5 * sqs
      rqs_over_sqs = rqs / sqs

      # Update beta
      SigmaBetaq = solve(rqs_over_sqs * SigmaBeta0_inv + WtW)
      muBetaq = SigmaBetaq %*% (SigmaBeta0_inv_muBeta0 + crossprod(W, muYStarq - Z %*% muuq - varphi %*% muThetaq))

      # Update muYStar
      muYStar = W %*% muBetaq + Z %*% muuq + varphi %*% muThetaq
      muYStarq = muYStar + dnorm(muYStar) / (pnorm(muYStar)^y * (pnorm(muYStar) - 1)^(1-y))

      lbnew = LB(y, W, Z, varphi, SigmaBetaq, Sigmauq, SigmaThetaq, muYStar, Sigmau0, rqs, sqs, rqt, sqt, muPsi, sigmaPsi2, r0t, s0t, r0s, s0s, , SigmaBeta0, muBetaq, muBeta0, w0, muThetaq)
      diff = lbnew - lbold
    }

    dif = (lbnew - lbold)/abs(lbnew)
    lbold = lbnew
    lbrecord = c(lbrecord, lbnew)
    cat("count: ", count, ", lbnew: ", lbnew, ", dif: ", dif, "\n")
  }
  list(muThetaq = muThetaq, SigmaThetaq = SigmaThetaq, muBetaq = muBetaq, SigmaBetaq = SigmaBetaq, rqt = rqt, sqt = sqt, rqs = rqs, sqs = sqs, muYStarq = muYStarq, muYStar = muYStar, lbrecord = lbrecord)
}


