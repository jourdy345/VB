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
  (-exp(0.5 * sigma2 * j^2 + mu * j) * dnorm(-mu/sigma - sigma * j) * ((mu/(2 * sigma^3)) - (j / (2 * sigma)))) + (0.5 * j^2 * exp(0.5 * sigma2 * j^2 + mu * j) * (1 - pnorm(-mu/sigma - (sigma * j)))) + (exp((0.5 * sigma2 * j^2) - (mu * j)) * dnorm(mu/sigma - (sigma * j)) * ((mu/(2 * sigma^3)) + (j/(2 * sigma)))) + (0.5 * j^2 * exp((0.5 * sigma2 * j^2) - (mu * j)) * (1 - pnorm(mu/sigma - sigma * j)))
}

der_S1_mu = function(sigma2, sigma, mu, w0) {
  -w0 * ( (-mu/sigma * sqrt(2/pi) * exp(-mu^2 / (2 * sigma2))) + (1 - 2 * pnorm(-mu/sigma)) + (2*mu/sigma * dnorm(-mu/sigma)) )
}

der_S1_sigma2 = function(sigma2, sigma, mu, w0) {
  -w0 * (((1/(2 * sigma)) + (0.5 * mu^2 / sigma^3) * sqrt(2/pi) * exp(-mu^2 / (2 * sigma2))) - (mu^2 / sigma^3 * dnorm(-mu/sigma)))
}

der_S2_mu = function(sigma2, sigma, sigtq, mutq, mu, rsoverssq, rtoverstq, w0, J) {
  ((-0.25 * J * (J + 1) / w0) * der_S1_mu(sigma2, sigma, mu, w0)) - (0.5 * rsoverssq * rtoverstq * sum(((diag(sigtq) + mutq^2) * der_Qj_mu(sigma2, sigma, mu, (1:J)))))
}

der_S2_sigma2 = function(sigma2, sigma, sigtq, mutq, mu, rsoverssq, rtoverstq, w0, J) {
  ((-0.25 * J * (J + 1) / w0) * der_S1_sigma2(sigma2, sigma, mu, w0)) - (0.5 * rsoverssq * rtoverstq * sum(((diag(sigtq) + mutq^2) * der_Qj_sigma2(sigma2, sigma, mu, (1:J)))))
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
  -0.5 * (tr(WtW %*% sigbq) + tr(varphitvarphi %*% sigtq)) + sum(log(((pnorm(tmp))^y) * ((1 - pnorm(tmp))^(1-y)))) - 0.5 * J * (log(2 * pi) - (digamma(rsq_half) - log(ssq_half) + digamma(rtq_half) - log(stq_half))) + (0.25 * J * (J + 1) * ((sigpsiq * sqrt(2 / pi) * exp(-mupsiq^2/(2*sigpsiq2))) + (mupsiq * (1 - 2 * pnorm(-mupsiq / sigpsiq))))) - (0.5 * rsoverssq * rtoverstq * sum((diag(sigtq) + mutq^2) * (exp( (0.5 * sigpsiq2 * ((1:J)^2)) + (mupsiq * (1:J)) ) * (1 - pnorm( - mupsiq / sigpsiq - (sigpsiq * (1:J)))) + exp((0.5 * sigpsiq2 * ((1:J)^2)) - (mupsiq * (1:J))) * (1 - pnorm(mupsiq/sigpsiq - (sigpsiq * (1:J))))))) + (0.5 * J * (1 + log(2 * pi) + determinant(sigtq)$modulus[1])) + (rt0_half * log(st0_half)) - lgamma(rt0_half) + ((rt0_half + 1) * (digamma(rtq_half) - log(stq_half))) - st0_half * rtoverstq + rtq_half + log(stq_half) + lgamma(rtq_half) - ((1 + rtq_half) * digamma(rtq_half)) + rs0_half * log(ssq_half) - lgamma(rs0_half) + ((rs0_half + 1) * (digamma(rsq_half) - log(ssq_half))) - ss0_half * rsoverssq + rsq_half + log(ssq_half) + lgamma(rsq_half) - ((1 + rsq_half) * digamma(rsq_half)) + 0.5 * ( p + digamma(rsq_half) - log(ssq_half)  + determinant(sigb0_inv_sigbq)$modulus[1] - (rsoverssq * (tr(sigb0_inv_sigbq) + sum((mubq - mub0) * (solve(sigb0, mubq - mub0)))))) + log(0.5 * w0) - w0 * ((sigpsiq * sqrt(2 / pi) * exp(-mupsiq^2/(2 * sigpsiq2))) + (mupsiq * (1 - 2 * pnorm(-mupsiq/sigpsiq)))) + 0.5 * (log(2*pi*sigpsiq2) - 1)
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
  sigpsiq2 = mupsi0^2 / 100
  sigpsiq = sqrt(sigpsiq2)
  mutq = rep(0, J)
  muystar = rep(0, n)
  sigbq = sigb0
  sigtq = diag(1, J)
  count = 0
  lbold = LB(y, x, W, sigbq, varphi, sigtq, mubq, mutq, sigpsiq2, mupsiq, rsq, ssq, rtq, stq, rs0, ss0, rt0, st0, sigb0, mub0, w0)
  lbnew = 0
  lbrecord = c(lbold)
  while (dif > tol) {
    count = count + 1
    # Update theta
    cat('sigpsiq: ', sigpsiq, '\n')
    cat('mupsiq: ', mupsiq, '\n')
    cat('rsoverssq: ', rsoverssq, '\n')
    cat('rtoverstq: ', rtoverstq, '\n')
    cat('MGF_foldedNormal: ', MGF_foldedNormal(sigpsiq2, sigpsiq, mupsiq, (1:J)), '\n')
    sigtq = solve(varphitvarphi + (rsoverssq * rtoverstq * diag(MGF_foldedNormal(sigpsiq2, sigpsiq, mupsiq, (1:J)))))
    mutq = sigtq %*% crossprod(varphi, muystar - W %*% mubq)

    # Update tau^2
    stq = st0 + (rtoverstq * (tr((sigtq + tcrossprod(mutq)) %*% diag(MGF_foldedNormal(sigpsiq2, sigpsiq, mupsiq, (1:J))))))
    stq_half = stq / 2
    rtoverstq = rtq / stq

    # Update sigma^2
    sigb0_inv_sigbq = solve(sigb0, sigbq)
    ssq = ss0 + (rtoverstq * (tr((sigtq + tcrossprod(mutq)) %*% diag(MGF_foldedNormal(sigpsiq2, sigpsiq, mupsiq, (1:J)))))) + tr(sigb0_inv_sigbq) + sum((mubq - mub0) * solve(sigb0, (mubq-mub0)))
    rsoverssq = rsq / ssq

    # Update beta
    sigbq = solve(WtW + rsoverssq * sigb0_inv)
    mubq = sigbq %*% (rsoverssq * sigb0_inv_sigbq + crossprod(W, muystar - varphi %*% mutq))

    # Update muystar
    muystar_tmp = W %*% mubq + varphi %*% mutq
    muystar = muystar_tmp + dnorm(muystar_tmp)/(((pnorm(muystar_tmp))^y)*((pnorm(muystar_tmp) - 1)^(1-y)))

    a = 1
    # Update psi
    sigpsiq2_old = sigpsiq2
    mupsiq_old = mupsiq
    sigpsiq2 = -0.5 / (der_S1_sigma2(sigpsiq2, sigpsiq, mupsiq, w0) + der_S2_sigma2(sigpsiq2, sigpsiq, sigtq, mutq, mupsiq, rsoverssq, rtoverstq, w0, J))
    # temp2 = sigpsiq2_new
    # mupsiq_old = mupsiq

    # print(der_S2_mu(sigpsiq2, sigpsiq, sigtq, mutq, mupsiq, rsoverssq, rtoverstq, w0, J) == der_S1_mu(sigpsiq2, sigpsiq, mupsiq, w0))
    # temp = der_S1_mu(sigpsiq2, sigpsiq, mupsiq, w0) + der_S2_mu(sigpsiq2, sigpsiq, sigtq, mutq, mupsiq, rsoverssq, rtoverstq, w0, J)
    # while (sigpsiq2_new < 0) {
    #   a = 2/3 * a
    #   sigpsiq2_new = 1/(1/sigpsiq2_old + a * (1/sigpsiq2_new - 1/sigpsiq2_old))
    # }
    # sigpsiq2 = sigpsiq2_new
    sigpsiq = sqrt(sigpsiq2)
    # mupsiq = mupsiq + a * sigpsiq2 * temp
    cat('der_S1_mu: ', der_S1_mu(sigpsiq2, sigpsiq, mupsiq, w0), '\n')
    cat('der_S2_mu: ', der_S2_mu(sigpsiq2, sigpsiq, sigtq, mutq, mupsiq, rsoverssq, rtoverstq, w0, J), '\n')
    mupsiq = mupsiq + a * sigpsiq2 * (der_S1_mu(sigpsiq2, sigpsiq, mupsiq, w0) + der_S2_mu(sigpsiq2, sigpsiq, sigtq, mutq, mupsiq, rsoverssq, rtoverstq, w0, J))

    lbnew = LB(y, x, W, sigbq, varphi, sigtq, mubq, mutq, sigpsiq2, mupsiq, rsq, ssq, rtq, stq, rs0, ss0, rt0, st0, sigb0, mub0, w0)
    lbfull <- lbnew
    dif = lbnew - lbold
    if (dif < 0) {
      a = 1
      dif_try <- dif
      while (dif_try < 0) {
        a = 0.5 * a
        sigpsiq2_try <- 1/(1/sigpsiq2_old + a*(1/sigpsiq2 - 1/sigpsiq2_old))
        sigpsiq_try <- sqrt(sigpsiq2_try)
        mupsiq_try <- sigpsiq2_try*(mupsiq_old/sigpsiq2_old + a * (mupsiq/sigpsiq2 - mupsiq_old/sigpsiq2_old ))
        lbnew <- LB(y, x, W, sigbq, varphi, sigtq, mubq, mutq, sigpsiq2, mupsiq, rsq, ssq, rtq, stq, rs0, ss0, rt0, st0, sigb0, mub0, w0)
        dif_try <- lbnew - lbfull
        # sigpsiq2_new = 1/(1/sigpsiq2 + a * (1/temp2 - 1/sigpsiq2))
        # sigpsiq2 = sigpsiq2_new
        # sigpsiq = sqrt(sigpsiq2)
        # mupsiq = mupsiq_old + a * sigpsiq2 * temp
      }

      sigpsiq2 <- sigpsiq2_try
      sigpsiq <- sigpsiq_try
      mupsiq <- mupsiq_try
      # while (sigpsiq2_new < 0) {
      #   a = 2/3 * a
      #   sigpsiq2_new = 1/(1/sigpsiq2 + a * (1/temp2 - 1/sigpsiq2))
      #   # sigpsiq2_try = 1/(1/sigpsiq2_old + step * (1/sigpsiq2 - 1/sigpsiq2_old))
      #   # sigpsiq_try = sqrt(sigpsiq2_try)
      #   # mupsiq_try = sigpsiq2_try * (mupsiq_old/sigpsiq2_old + step*(mupsiq/sigpsiq2 - mupsiq_old/sigpsiq2_old))
      #   # lbnew = LB(y, x, W, sigbq, varphi, sigtq, mubq, mutq, sigpsiq2_try, mupsiq_try, rsq, ssq, rtq, stq, rs0, ss0, rt0, st0, sigb0, mub0, w0)
      #   # dif_try = lbnew - lbtmp
      #   # cat('sigpsiq2_try: ', sigpsiq2_try, ', mupsiq_try: ', mupsiq_try, '\n')
      # }
      # sigpsiq2 = sigpsiq2_new
      # sigpsiq = sqrt(sigpsiq2)
      # mupsiq = mupsiq_old + a * sigpsiq2 * temp

      # Update theta
      # sigtq = solve(varphitvarphi + (rsoverssq * rtoverstq * diag(MGF_foldedNormal(sigpsiq2, sigpsiq, mupsiq, (1:J)))))
      # mutq = sigtq %*% crossprod(varphi, muystar - W %*% mubq)

      # # Update tau^2
      # stq = st0 + (rtoverstq * (tr((sigtq + tcrossprod(mutq)) %*% diag(MGF_foldedNormal(sigpsiq2, sigpsiq, mupsiq, (1:J))))))
      # stq_half = stq / 2
      # rtoverstq = rtq / stq

      # # Update sigma^2
      # sigb0_inv_sigbq = solve(sigb0, sigbq)
      # ssq = ss0 + (rtoverstq * (tr((sigtq + tcrossprod(mutq)) %*% diag(MGF_foldedNormal(sigpsiq2, sigpsiq, mupsiq, (1:J)))))) + tr(sigb0_inv_sigbq) + sum((mubq - mub0) * solve(sigb0, (mubq-mub0)))
      # rsoverssq = rsq / ssq

      # # Update beta
      # sigbq = solve(rsoverssq * (WtW + sigb0_inv))
      # mubq = sigbq %*% (rsoverssq * sigb0_inv_sigbq + crossprod(W, muystar - varphi %*% mutq))

      # # Update muystar
      # muystar_tmp = W %*% mubq + varphi %*% mutq
      # muystar = muystar_tmp + dnorm(muystar_tmp)/(((pnorm(muystar_tmp))^y)*((pnorm(muystar_tmp) - 1)^(1-y)))

      # lbnew = LB(y, x, W, sigbq, varphi, sigtq, mubq, mutq, sigpsiq2, mupsiq, rsq, ssq, rtq, stq, rs0, ss0, rt0, st0, sigb0, mub0, w0)
      # diff = lbnew - lbold
    }

    dif = (lbnew - lbold)/abs(lbnew)
    lbold = lbnew
    lbrecord = c(lbrecord, lbnew)
    cat("count: ", count, ", lbnew: ", lbnew, ", dif: ", dif, "\n")
  }
  list(mutq = mutq, sigtq = sigtq, mubq = mubq, sigbq = sigbq, rtq = rtq, stq = stq, rsq = rsq, ssq = ssq, muystar = muystar, lbrecord = lbrecord)
}

sim_cosine_probit = function(FUN, J, draw = TRUE) {


  ##########################################################

  ################    Start simulation    ##################

  ##########################################################

  N_ = 1000
  x = runif(N_, 0, 1)
  W = rep(1, times = N_)
  mub0 = 0
  y = rbinom(length(x),1,pnorm(FUN(x) + W))
  sigb0 = matrix(1,nrow=1,ncol=1)
  fit = VB(y, x, W, mub0, sigb0, w0 = 1, rs0 = 0.01, ss0 = 0.01, rt0 = 0.01, st0 = 0.01, mupsi0 = 0.01, J = J, tol = 1.0e-06)
  # categorize = function(x) {
  #   y = 0
  #   if (x < 0) {
  #     y = 0
  #   } else {
  #     y = 1
  #   }
  #   y
  # }

  varphi = sqrt(2)*cos(outer(x,pi*(1:J)))
  unknown_f = varphi %*% fit$mutq + W %*% fit$mubq
  # res = sapply(fit$muystar, categorize)
  if (draw == TRUE) {
    # plot(y-mean(y) ~ xobs, xlab = 'index', ylab = 'observed/fitted', main = 'Simulation result', type = 'p')
    ord = order(x)
    plot(x[ord], unknown_f[ord], col = 'purple', xlab = 'index', ylab = 'observed/fitted', main = 'Simulation result', ylim = c(-10, 10), type = 'l')
    # lines(xobs[ord], fixed[ord], col = 'purple')
    curve(FUN, from = 0, to = 1, col = 'red', lty = 2, add = TRUE)
    legend("topright", legend = c('fitted values', 'true function'), col = c('purple', 'red'), lty = c(1, 2), bg = 'gray95')
    return(list(fit = fit, y = y, x = x))
  } else {
    return(list(fit = fit, y = y, x = x, W = W, unknown_f = unknown_f))
  }
}

# sim_cosine_probit = function(FUN, J) {
#   der_Qj_mu = function(sigma2, sigma, mu, j) {
#     (exp((0.5 * sigma2 * j^2) + (mu * j)) * dnorm(-mu/sigma - sigma * j) / sigma) + (j * exp((0.5 * sigma2 * j^2) + (mu * j)) * (1 - pnorm((-mu/sigma) - (sigma * j)))) - (exp(0.5 * sigma2 * j^2 - mu * j) * dnorm(mu/sigma - sigma * j) / sigma) - (j * exp((0.5 * sigma2 * j^2) - (mu * j)) * (1 - pnorm(mu/sigma - sigma * j)))
#   }
#   der_Qj_sigma2 = function(sigma2, sigma, mu, j) {
#     (-exp(0.5 * sigma2 * j^2 + mu * j) * dnorm(-mu/sigma - sigma * j) * ((mu/(2 * sigma^3)) - (j / (2 * sigma2)))) + (0.5 * j^2 * exp(0.5 * sigma2 * j^2 + mu * j) * (1 - pnorm(-mu/sigma - (sigma * j)))) + (exp((0.5 * sigma2 * j^2) - (mu * j)) * dnorm(mu/sigma - (sigma * j)) * ((mu/(2 * sigma^3)) + (j/(2 * sigma2)))) + (0.5 * j^2 * exp((0.5 * sigma2 * j^2) - (mu * j)) * (1 - pnorm(mu/sigma - sigma * j)))
#   }

#   der_S1_mu = function(sigma2, sigma, mu, w0) {
#     -w0 * ( (-mu/sigma * sqrt(2/pi) * exp(-mu^2 / (2 * sigma2))) + (1 - 2 * pnorm(-mu/sigma)) + (2*mu/sigma * dnorm(-mu/sigma)) )
#   }

#   der_S1_sigma2 = function(sigma2, sigma, mu, w0) {
#     -w0 * (((1/(2 * sigma)) + (0.5 * mu^2 / sigma^3) * sqrt(2/pi) * exp(-mu^2 / (2 * sigma2))) - (mu^2 / sigma^3 * dnorm(-mu/sigma)))
#   }

#   der_S2_mu = function(sigma2, sigma, sigtq, mutq, mu, rsoverssq, rtoverstq, w0, J) {
#     ((-0.25 * J * (J + 1) / w0) * der_S1_mu(sigma2, sigma, mu, w0)) - (0.5 * rsoverssq * rtoverstq * sum(((diag(sigtq) + mutq^2) * der_Qj_mu(sigma2, sigma, mu, (1:J)))))
#   }

#   der_S2_sigma2 = function(sigma2, sigma, sigtq, mutq, mu, rsoverssq, rtoverstq, w0, J) {
#     ((-0.25 * J * (J + 1) / w0) * der_S1_sigma2(sigma2, sigma, mu, w0)) - (0.5 * rsoverssq * rtoverstq * sum(((diag(sigtq) + mutq^2) * der_Qj_sigma2(sigma2, sigma, mu, (1:J)))))
#   }
#   Qj = function(mu, sigma2, j) {
#     term1 = 0.5 * sigma2 * j^2
#     term2 = mu * j
#     term3 = mu/sqrt(sigma2)
#     term4 = sqrt(sigma2) * j
#     (exp(term1 + term2) * (1 - pnorm(-(term3 + term4)))) + (exp(term1 - term2) * (1 - pnorm(term3 - term4)))
#   }
#   S1 = function(mu, sigma2, w0) -w0 * ((sqrt(sigma2 * 2 / pi) * exp(-mu^2 / (2 * sigma2))) + (mu * (1 - 2 * pnorm(-mu/sqrt(sigma2)))))
#   S2 = function(mupsiq, sigpsiq2, w0, J, rtq, stq, rsq, ssq, sigtq, mutq) ((-0.5 * rtq/stq * rsq/ssq) * sum((diag(sigtq) + mutq^2) * Qj(mupsiq, sigpsiq2, 1:J))) - (J*(J+1)/(4*w0) * S1(mupsiq, sigpsiq2, w0))

#   LB = function(y, W, varphi, mubq, mutq, sigb0, sigbq, rsq, ssq, rs0, ss0, rtq, stq, rt0, st0, mub0, mupsiq, sigpsiq2, w0, J) {
#     p = dim(W)[2]
#     st0_half = st0/2
#     ss0_half = ss0/2
#     rt0_half = rt0/2
#     rs0_half = rs0/2
#     stq_half = stq/2
#     ssq_half = ssq/2
#     rtq_half = rtq/2
#     rsq_half = rsq/2
#     tau_ratio0 = rt0/st0
#     sig_ratio0 = rs0/ss0
#     tau_ratioq = rtq/stq
#     sig_ratioq = rsq/ssq
#     WtW = crossprod(W)
#     PtP = crossprod(varphi)
#     t1 = W %*% mubq + varphi %*% mutq
#     t2 = solve(sigb0, sigbq)
#     (-0.5 * sum(WtW * sigbq) + sum(PtP * sigtq)) + sum(log((pnorm(t1)^y) * ((1-pnorm(t1))^(1-y))))
#     + 0.5 * (J * (digamma(rsq_half) - log(ssq_half) + digamma(rtq_half) - log(stq_half) + 1) + determinant(sigtq, logarithm = TRUE)$modulus[1] + determinant(t2, logarithm = TRUE)$modulus[1])
#     + rt0_half * log(st0_half) - lgamma(rt0_half) + (rt0_half + 1) * (digamma(rtq_half) - log(stq_half)) - st0_half * tau_ratioq + rtq_half + log(stq_half) + lgamma(rtq_half)
#     - (1 + rtq_half) * digamma(rtq_half)
#     +rs0_half * log(ss0_half) - lgamma(rs0_half) + (rs0_half + 1) * (digamma(rsq_half)-log(ssq_half)) - ss0_half * sig_ratioq + rsq_half +log(ssq_half) + lgamma(rsq_half)
#     - (1 + rsq_half) * digamma(rsq_half)
#     + 0.5 * (p + digamma(rsq_half) - log(ssq_half) - sig_ratioq * (sum(diag(t2)) + sum((mubq - mub0) * solve(sigb0, mubq - mub0))) + log(2*pi*sigpsiq2) + 1) + log(w0/2)
#     + S1(mupsiq, sigpsiq2, w0) + S2(mupsiq, sigpsiq2, w0, J, rtq, stq, rsq, ssq, sigtq, mutq)
#   }

#   VB = function(y, W, x, st0, st0, rs0, ss0, mub0, sigb0, w0, J, tol = 1.0e-6) {
#     WtW = crossprod(W)
#     PtP = crossprod(varphi)
#     p = dim(W)[2]
#     # Initialize variational parameters
#     rtq = rt0 + J
#     stq = 0.01
#     tau_ratioq = rtq/stq
#     rsq = rs0 + J + p
#     ssq = 0.01
#     sig_ratioq = rsq/ssq
#     mutq = rep(0, J)
#     sigtq = diag(1, J)
#     sigbq = diag(1, p)
#     mubq = rep(0, p)
#     sigpsiq2 = 0.01
#     mupsiq = 0.1
#     muystar = rep(0, length(y))

#     sigb0_inv = solve(sigb0)

#     dif = tol + 1
#     lbold = 10e-7
#     lbnew = 0
#     lbrecord = c()
#     count = 0
#     while (dif > tol) {
#       count = count + 1
#       a = 1
#       # Update psi
#       sigpsiq2_old = sigpsiq2
#       sigpsiq2_new = -0.5 / ((der_S1_sigma2(sigpsiq2, sigpsiq, mupsiq, w0) + der_S2_sigma2(sigpsiq2, sigpsiq, sigtq, mutq, mupsiq, sig_ratioq, tau_ratioq, w0, J)))
#       temp2 = sigpsiq2_new
#       mupsiq_old = mupsiq
#       temp = der_S1_mu(sigpsiq2, sigpsiq, mupsiq, w0) + der_S2_mu(sigpsiq2, sigpsiq, sigtq, mutq, mupsiq, sig_ratioq, tau_ratioq, w0, J)
#       while (sigpsi2_new < 0) {
#         a = 2/3 * a
#         sigpsiq2_new = 1/(1/sigpsiq2_old + a * (1/sigpsiq2_new - 1/sigpsiq2_old))
#       }
#       sigpsiq2 = sigpsiq2_new
#       sigpsiq = sqrt(sigpsiq2)
#       mupsiq = mupsiq + a * sigpsiq2 * temp

#       # Update theta
#       sigtq = solve(PtP + sig_ratioq * tau_ratioq * diag(Qj(mupsiq, sigpsiq2, 1:J)))
#       mutq = sigtq %*% crossprod(varphi, muystar - W %*% mubq)
      
#       # Update tau
#       stq = st0 + sig_ratioq * sum((diag(sigtq) + mutq^2) * Qj(mupsiq, sigpsiq2, 1:J))
#       tau_ratioq = rtq/stq

#       # Update sigma
#       ssq = ss0 + tau_ratioq * sum((diag(sigtq) + mutq^2) * Qj(mupsiq, sigpsiq2, 1:J)) + sum(sigb0_inv * sigbq) + sum((mubq - mub0) * (sigb0_inv %*% (mubq - mub0)))
#       sig_ratioq = rsq/ssq
      
#       # Update beta
#       sigbq = solve(sig_ratioq * (Wtw + sigb0_inv))
#       mubq = sig_ratioq * sigbq %*% (sigb0_inv %*% mub0 + crossprod(W, muystar - varphi %*% mutq))

#       # Update muystar

#     }


#   }
# }