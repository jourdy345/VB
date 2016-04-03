tr = function(X) sum(diag(X))
E_foldedNormal = function(mu, sigma2, sigma) {
  (sigma * sqrt(2 / pi) * exp(-mu^2/(2 * sigma2))) + (mu * (1 - 2 * pnorm(-mu / sigma)))
}
MGF_foldedNormal = function(sigma2, sigma, mu, t) {
  (exp( (0.5 * sigma2 * t^2) + (mu * t)) * (1 - pnorm(-mu/sigma - (sigma * t)))) + (exp((0.5 * sigma2 * t^2) - (mu * t)) * (1 - pnorm(mu/sigma - (sigma * t))))
}

der_Qj_mu = function(sigma2, sigma, mu, j) {
  (exp((0.5 * sigma2 * j^2) + (mu * j)) * dnorm(-mu/sigma - sigma * j) / sigma) + (j * exp((0.5 * sigma2 * j^2) + (mu * j)) * (1 - pnorm((-mu/sigma) - (sigma * j)))) - (exp(0.5 * sigma2 * j^2 - mu * j) * dnorm(mu/sigma - sigma * j) / sigma) - (j * exp((0.5 * sigma2 * j^2) - (mu * j)) * (1 - pnorm(mu/sigma - sigma * j)))
}
der_Qj_sigma2 = function(sigma2, sigma, mu, j) {
  (-exp(0.5 * sigma2 * j^2 + mu * j) * dnorm(-mu/sigma - sigma * j) * ((mu/(2 * sigma^3)) - (j / (2 * sigma)))) + (0.5 * j^2 * exp(0.5 * sigma2 * j^2 + mu * j) * (1 - pnorm(-mu/sigma - (sigma * j)))) + (exp((0.5 * sigma2 * j^2) - (mu * j)) * dnorm(mu/sigma - (sigma * j)) * ((mu/(2 * sigma^2)) + (j/(2 * sigma)))) + (0.5 * j^2 * exp((0.5 * sigma2 * j^2) - (mu * j)) * (1 - pnorm(mu/sigma - sigma * j)))
}

der_S1_mu = function(sigma2, sigma, mu, w0) {
  -w0 * ( (-mu/sigma * sqrt(2/pi) * exp(-mu^2 / (2 * sigma2))) + (1 - 2 * pnorm(-mu/sigma)) + (2*mu/sigma * dnorm(-mu/sigma)) )
}

der_S1_sigma2 = function(sigma2, sigma, mu, w0) {
  -w0 * (((1/(2 * sigma)) + (0.5 * mu^2 / sigma^3) * sqrt(2/pi) * exp(mu^2 / (2 * sigma2))) - (mu^2 / sigma^3 * dnorm(-mu/sigma)))
}

der_S2_mu = function(sigma2, sigma, sigtq, mutq, mu, rsoverssq, rtoverstq, w0, J) {
  ((-0.25 * J * (J + 1) / w0) * der_S1_mu(sigma2, sigma, mu, w0)) - (0.5 * rsoverssq * rtoverstq * sum((diag(sigtq) + mutq^2) * der_Qj_mu(sigma2, sigma, mu, (1:J))))
}

der_S2_sigma2 = function(sigma2, sigma, sigtq, mutq, mu, rsoverssq, rtoverstq, w0, J) {
  ((-0.25 * J * (J + 1) / w0) * der_S1_sigma2(sigma2, sigma, mu, w0)) - (0.5 * rsoverssq * rtoverstq * sum((diag(sigtq) + mutq^2) * der_Qj_sigma2(sigma2, sigma, mu, (1:J))))
}


LB = function(y, x, W, sigbq, varphi, sigtq, mubq, mutq, sigpsiq2, mupsiq, rsq, ssq, rtq, stq, rs0, ss0, rt0, st0, sigb0, mub0, w0) {
  J = length(mutq)
  p = length(mubq)
  WtW = crossprod(W)
  varphitvarphi = crossprod(varphi)
  tmp = W %*% mubq + varphi %*% mutq
  sigpsiq = sqrt(sigpsiq2)
  rt0_half = 0.5 * rt0
  st0_half = 0.5 * st0
  rtq_half = 0.5 * rtq
  stq_half = 0.5 * stq
  rtoverstq = rtq / stq
  rs0_half = 0.5 * rs0
  ss0_half = 0.5 * ss0
  rsq_half = 0.5 * rsq
  ssq_half = 0.5 * ssq
  rsoverssq = rsq / ssq
  sigb0_inv_sigbq = solve(sigb0, sigbq)
  -0.5 * (tr(WtW %*% sigbq) + tr(varphitvarphi %*% sigtq)) + sum(log(((pnorm(tmp))^y) * ((1 - pnorm(tmp))^(1-y)))) -0.5 * J * (log(2 * pi) - (digamma(rsq_half) - log(ssq_half) + digamma(rtq_half) - log(stq_half))) + 0.25 * J * (J + 1) * (sigpsiq * sqrt(2 / pi) * exp(-mupsiq^2/(2*sigpsiq2)) + mupsiq * (1 - 2 * pnorm(-mupsiq / sigpsiq))) - 0.5 * rsoverssq * rtoverstq * sum((diag(sigtq) + mutq^2) * (exp( (0.5 * sigpsiq2 * (1:J)) + (mupsiq * (1:J)) ) * (1 - pnorm( - mupsiq / sigpsiq - (sigpsiq * (1:J)))) + exp((0.5 * sigpsiq2 * (1:J)) - (mupsiq * (1:J))) * (1 - pnorm(mupsiq/sigpsiq - (sigpsiq * (1:J)))))) + (rt0_half * log(st0_half)) - lgamma(rt0_half) + ((rt0_half + 1) * (digamma(rtq_half) - log(stq_half))) - st0_half * rtoverstq + rtq_half + log(stq_half) + lgamma(rtq_half) - ((1 + rtq_half) * digamma(rtq_half)) + rs0_half * log(ssq_half) - lgamma(rs0_half) + ((rs0_half + 1) * (digamma(rsq_half) - log(ssq_half))) - ss0_half * rsoverssq + rsq_half + log(ssq_half) + lgamma(rsq_half) - ((1 + rsq_half) * digamma(rsq_half)) + 0.5 * ( p + digamma(rsq_half) - log(ssq_half)  + determinant(sigb0_inv_sigbq)$modulus[1] - (rsoverssq * (tr(sigb0_inv_sigbq) + sum((mubq - mub0) * (solve(sigb0, mubq - mub0)))))) + log(0.5 * w0) - w0 * ((sigpsiq * sqrt(2 / pi) * exp(-mupsiq^2/(2 * sigpsiq2))) + (mupsiq * (1 - 2 * pnorm(-mupsiq/sigpsiq)))) + 0.5 * (log(2*pi*sigpsiq2) - 1)
}

VB = function(y, x, W, mub0, sigb0, w0, rs0, ss0, rt0, st0, mupsi0, J, tol = 1.0e-06) {
  n = length(y)
  if (!is.matrix(W)) W = as.matrix(W)
  p = dim(W)[2]
  dif = tol + 1
  varphi = sqrt(2)*cos(outer(x,pi*(1:J)))
  varphitvarphi = crossprod(varphi)
  WtW = crossprod(W)
  sigb0_inv = solve(sigb0)
  sigb0_inv_mub0 = solve(sigb0, mub0)
  # Initialize variational parameters
  rsq = rs0 + J + p
  ssq = ss0
  rtq = rt0 + J
  stq = st0
  st0_half = st0/2
  ss0_half = ss0/2
  rt0_half = rt0/2
  rs0_half = rs0/2
  rtoverstq = rtq / stq
  rsoverssq = rsq / ssq
  mubq = mub0
  mupsiq = mupsi0
  sigpsiq2 = mupsi0^2 / 10
  sigpsiq = sqrt(sigpsiq2)
  mutq = rnorm(J, 0, st0_half/(rt0_half - 1) * ss0_half/(rs0_half - 1) * exp(-(1:J) * abs(mupsi0)))
  muystar = W %*% mubq + varphi %*% mutq
  sigbq = sigb0
  sigtq = solve((rsoverssq * varphitvarphi) + (rsoverssq * rtoverstq * diag(MGF_foldedNormal(sigpsiq2, sigpsiq, mupsiq, (1:J)))))

  count = 0
  lbold = LB(y, x, W, sigbq, varphi, sigtq, mubq, mutq, sigpsiq2, mupsiq, rsq, ssq, rtq, stq, rs0, ss0, rt0, st0, sigb0, mub0, w0)
  lbrecord = c(lbold)
  while (dif > tol) {
    count = count + 1
    # Update theta
    sigtq = solve((rsoverssq * varphitvarphi) + (rsoverssq * rtoverstq * diag(MGF_foldedNormal(sigpsiq2, sigpsiq, mupsiq, (1:J)))))
    mutq = sigtq %*% crossprod(varphi, muystar - W %*% mubq)

    # Update tau^2
    stq = st0 + (rtoverstq * (tr((sigtq + tcrossprod(mutq)) %*% diag(MGF_foldedNormal(sigpsiq2, sigpsiq, mupsiq, (1:J))))))
    stq_half = stq / 2
    rtoverstq = rtq / stq

    # Update sigma^2
    sigb0_inv_sigbq = solve(sigb0, sigbq)
    ssq = ss0 + (rtoverstq * (tr((sigtq + tcrossprod(mutq)) %*% diag(MGF_foldedNormal(sigpsiq2, sigpsiq, mupsiq, (1:J)))))) + tr(sigb0_inv_sigbq) + sum((mubq - mub0) * solve(sigb0, (mubq-mub0)))

    # Update beta
    sigbq = solve(rsoverssq * (WtW + sigb0_inv))
    mubq = rsoverssq * sigbq %*% (sigb0_inv_sigbq + crossprod(W, muystar - varphi %*% mutq))
    lbtmp = LB(y, x, W, sigbq, varphi, sigtq, mubq, mutq, sigpsiq2, mupsiq, rsq, ssq, rtq, stq, rs0, ss0, rt0, st0, sigb0, mub0, w0)
    
    # Reserve psi for step-halving
    sigpsi2_old = sigpsiq2
    mupsiq_old = mupsiq

    # Update psi
    sigpsiq2 = -0.5 / (der_S1_sigma2(sigpsiq2, sigpsiq, mupsiq, w0) + der_S2_sigma2(sigpsiq2, sigpsiq, sigtq, mutq, mupsiq, rsoverssq, rtoverstq, w0, J))
    mupsiq = mupsiq + sigpsiq2 * (der_S1_mu(sigma2, sigma, mu, w0) + der_S2_mu(sigpsiq2, sigpsiq, sigtq, mutq, mupsiq, rsoverssq, rtoverstq, w0, J))
    sigpsiq = sqrt(sigpsiq2)
    lbnew = LB(y, x, W, sigbq, varphi, sigtq, mubq, mutq, sigpsiq2, mupsiq, rsq, ssq, rtq, stq, rs0, ss0, rt0, st0, sigb0, mub0, w0)
    dif = lbnew - lbtmp
    if (dif < 0) {
      step = 1
      dif_try = dif
      while (dif_try < 0) {
        step = 0.5 * step
        sigpsiq2_try = 1/(1/sigpsiq2_old + step * (1/sigpsiq2 - 1/sigpsiq2_old))
        sigpsiq_try = sqrt(sigpsiq2_try)
        mupsiq_try = sigpsiq2_try * (mupsiq_old/sigpsiq2_old + step*(mupsiq/sigpsiq2 - mupsiq_old/sigpsiq2_old))
        lbnew = LB(y, x, W, sigbq, varphi, sigtq, mubq, mutq, sigpsiq2_try, mupsiq_try, rsq, ssq, rtq, stq, rs0, ss0, rt0, st0, sigb0, mub0, w0)
        dif_try = lbnew - lbtmp
      }
      sigpsiq2 = sigpsiq2_try
      sigpsiq = sigpsiq_try
      mupsiq = mupsiq_try
    }

    dif = (lbnew - lbold)/abs(lbnew)
    lbold = lbnew
    lbrecord = c(lbrecord, lbold)
    cat("count: ", count, ", lbnew: ", lbnew, ", dif: ", dif, "\n")
  }
  list(mutq = mutq, sigtq = sigtq, mubq = mubq, sigbq = sigbq, rtq = rtq, stq = stq, rsq = rsq, ssq = ssq, muystar = muystar, lbrecord = lbrecord)
}