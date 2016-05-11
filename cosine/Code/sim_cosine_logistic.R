sim_cosine_logistic = function(FUN, J) {
  lambda_xi = function(x) -tanh(x/2)/(4*x)
  Psi_xi = function(x) x/2 - log(1 + exp(x)) + x*tanh(x/2)/4

  der_Qj_mu = function(sigma2, sigma, mu, j) {
    (exp((0.5 * sigma2 * j^2) + (mu * j)) * dnorm(-mu/sigma - sigma * j) / sigma) + (j * exp((0.5 * sigma2 * j^2) + (mu * j)) * (1 - pnorm((-mu/sigma) - (sigma * j)))) - (exp(0.5 * sigma2 * j^2 - mu * j) * dnorm(mu/sigma - sigma * j) / sigma) - (j * exp((0.5 * sigma2 * j^2) - (mu * j)) * (1 - pnorm(mu/sigma - sigma * j)))
  }
  der_Qj_sigma2 = function(sigma2, sigma, mu, j) {
    (-exp(0.5 * sigma2 * j^2 + mu * j) * dnorm(-mu/sigma - sigma * j) * ((mu/(2 * sigma^3)) - (j / (2 * sigma2)))) + (0.5 * j^2 * exp(0.5 * sigma2 * j^2 + mu * j) * (1 - pnorm(-mu/sigma - (sigma * j)))) + (exp((0.5 * sigma2 * j^2) - (mu * j)) * dnorm(mu/sigma - (sigma * j)) * ((mu/(2 * sigma^3)) + (j/(2 * sigma2)))) + (0.5 * j^2 * exp((0.5 * sigma2 * j^2) - (mu * j)) * (1 - pnorm(mu/sigma - sigma * j)))
  }

  der_S1_mu = function(sigma2, sigma, mu, w0) {
    -w0 * ( (-mu/sigma * sqrt(2/pi) * exp(-mu^2 / (2 * sigma2))) + (1 - 2 * pnorm(-mu/sigma)) + (2*mu/sigma * dnorm(-mu/sigma)) )
  }

  der_S1_sigma2 = function(sigma2, sigma, mu, w0) {
    -w0 * (((1/(2 * sigma)) + (0.5 * mu^2 / sigma^3) * sqrt(2/pi) * exp(-mu^2 / (2 * sigma2))) - (mu^2 / sigma^3 * dnorm(-mu/sigma)))
  }

  der_S2_mu = function(sigma2, sigma, sigtq, mutq, mu, rtoverstq, w0, J) {
    ((-0.25 * J * (J + 1) / w0) * der_S1_mu(sigma2, sigma, mu, w0)) - (0.5 * rtoverstq * sum(((diag(sigtq) + mutq^2) * der_Qj_mu(sigma2, sigma, mu, (1:J)))))
  }

  der_S2_sigma2 = function(sigma2, sigma, sigtq, mutq, mu, rtoverstq, w0, J) {
    ((-0.25 * J * (J + 1) / w0) * der_S1_sigma2(sigma2, sigma, mu, w0)) - (0.5 * rtoverstq * sum(((diag(sigtq) + mutq^2) * der_Qj_sigma2(sigma2, sigma, mu, (1:J)))))
  }
  Qj = function(mu, sigma2, j) {
    term1 = 0.5 * sigma2 * j^2
    term2 = mu * j
    term3 = mu/sqrt(sigma2)
    term4 = sqrt(sigma2) * j
    (exp(term1 + term2) * (1 - pnorm(-(term3 + term4)))) + (exp(term1 - term2) * (1 - pnorm(term3 - term4)))
  }
  S1 = function(mu, sigma2, w0) -w0 * ((sqrt(sigma2 * 2 / pi) * exp(-mu^2 / (2 * sigma2))) + (mu * (1 - 2 * pnorm(-mu/sqrt(sigma2)))))
  S2 = function(mupsiq, sigpsiq2, w0, J, rtq, stq, sigtq, mutq) (-0.5 * rtq/stq * sum((diag(sigtq) + mutq^2) * Qj(mupsiq, sigpsiq2, 1:J))) - (J*(J+1)/(4*w0) * S1(mupsiq, sigpsiq2, w0))
  
  LB = function(mutq, varphi, xi, sigtq, rtq, stq, rt0, st0, w0, mupsiq, sigpsiq2) {
    J = dim(varphi)[2]
    t1 = as.vector(varphi %*% mutq)
    lambda = lambda_xi(xi)
    half_rtq = rtq / 2
    half_stq = stq / 2
    half_rt0 = rt0 / 2
    half_st0 = st0 / 2
    ratioq = rtq / stq
    # varphi = sqrt(2)*cos(outer(x,pi*(1:J)))
    crossprod(t1, y - 0.5) + sum(lambda * diag(varphi %*% tcrossprod(sigtq, varphi))) + sum(t1^2 * lambda) + sum(Psi_xi(xi))
    + 0.5 * J * (digamma(half_rtq) - log(half_stq) - log(2*pi)) + S1(mupsiq, sigpsiq2, w0) + S2(mupsiq, sigpsiq2, w0, J, rtq, stq, sigtq, mutq)
    + half_rt0 * log(half_st0) - lgamma(half_rt0) + half_rt0 * (digamma(half_rtq) - log(half_stq)) - half_st0 * ratioq + log(w0/2)
    + half_rtq + lgamma(half_rtq) - half_rtq * digamma(half_rtq) + 0.5 * (J * (1 + log(2*pi)) + determinant(sigtq)$modulus[1] + log(2*pi*sigpsiq2) + 1)
  }


  VB = function(y, x, rt0, st0, w0, J, tol = 1.0e-6) {
    varphi = sqrt(2)*cos(outer(x,pi*(1:J)))

    # Initialize variational parameters
    rtq = rt0 + J
    stq = rt0
    ratioq = rtq/stq
    sigtq = diag(1, J)
    mutq = rep(0, J)
    xi = rep(1, length(y))
    mupsiq = 1
    sigpsiq2 = 1/100
    sigpsiq = sqrt(sigpsiq2)

    lbold = -10e7
    lbtmp = 0
    dif = tol + 1
    lbrecord = c()
    count = 0
    while(dif > tol) {
      count = count + 1
      a = 1
      cat('count: ', count, '\n')
      # Update psi
      sigpsiq2_old = sigpsiq2
      sigpsiq2_new = -0.5 / (der_S1_sigma2(sigpsiq2, sigpsiq, mupsiq, w0) + der_S2_sigma2(sigpsiq2, sigpsiq, sigtq, mutq, mupsiq, ratioq, w0, J))
      temp2 = sigpsiq2_new
      mupsiq_old = mupsiq
      temp = der_S1_mu(sigpsiq2, sigpsiq, mupsiq, w0) + der_S2_mu(sigpsiq2, sigpsiq, sigtq, mutq, mupsiq, ratioq, w0, J)
      while (sigpsiq2_new < 0) {
        a = 2/3 * a
        sigpsiq2_new = 1/(1/sigpsiq2_old + a * (1/sigpsiq2_new - 1/sigpsiq2_old))
      }
      sigpsiq2 = sigpsiq2_new
      sigpsiq = sqrt(sigpsiq2)
      mupsiq = mupsiq + a * sigpsiq2 * temp
  
      # Update theta
      # cat('theta part1: ', crossprod(varphi, varphi * lambda_xi(xi)), '\n')
      sigtq = solve(crossprod(varphi, varphi * (-2 * lambda_xi(xi))) + ratioq * diag(Qj(mupsiq, sigpsiq2, 1:J)))
      mutq = sigtq %*% crossprod(varphi, y - 0.5)

      # Update tau
      stq = st0 + sum((diag(sigtq) + mutq^2) * Qj(mupsiq, sigpsiq2, 1:J))
      ratioq = rtq/stq
      # Update xi
      xi2 = diag(varphi %*% tcrossprod(sigtq + tcrossprod(mutq), varphi))
      xi = sqrt(xi2)

      # Check LB
      lbtmp = LB(mutq, varphi, xi, sigtq, rtq, stq, rt0, st0, w0, mupsiq, sigpsiq2)
      diff = lbtmp - lbold
      cat('lbtmp: ', lbtmp, '\n')
      a = 1
      if (diff < 0) {
        sigpsiq2_new = temp2
        
        while(sigpsiq2_new < 0) {
          a = 2/3 * a
          sigpsiq2_new = 1/(1/sigpsiq2 + a * (1/temp2 - 1/sigpsiq2))
        }

        sigpsiq2 = sigpsiq2_new
        sigpsiq = sqrt(sigpsiq2)
        mupsiq = mupsiq_old + a * sigpsiq2 * temp
        cat('mupsiq: ', mupsiq, '\n')
        
        # Update tau
        stq = st0 + sum((diag(sigtq) + mutq^2) * Qj(mupsiq, sigpsiq2, 1:J))
        ratioq = rtq/stq
        # Update theta
        # cat('lambda: ', lambda_xi(xi), '\n')
        # cat('\n part1: ', (crossprod(varphi, (varphi * (-2*lambda_xi(xi))))), '\n')
        # cat('\n last: ', Qj(mupsiq, sigpsiq2, 1:J), '\n')
        # cat('\n inside: ', eigen((crossprod(varphi, (varphi * (-2*lambda_xi(xi))))) + ratioq * diag(Qj(mupsiq, sigpsiq2, 1:J)))$values, '\n')
        sigtq = solve((crossprod(varphi, (varphi * (-2*lambda_xi(xi))))) + ratioq * diag(Qj(mupsiq, sigpsiq2, 1:J)))
        mutq = sigtq %*% crossprod(varphi,  y - 0.5)

        # Update xi
        xi = sqrt(diag(varphi %*% tcrossprod(sigtq + tcrossprod(mutq), varphi)))
        lbtmp = LB(mutq, varphi, xi, sigtq, rtq, stq, rt0, st0, w0, mupsiq, sigpsiq2)
        diff = lbtmp - lbold
        cat('latter diff: ', diff, '\n')
      }
      

      dif = (lbtmp - lbold)/abs(lbtmp)
      cat('dif: ', dif, '\n')
      lbold = lbtmp
      lbrecord = c(lbrecord, lbtmp)
    }
    list(xi = xi, sigtq = sigtq, mutq = mutq, rtq = rtq, stq = stq, lbrecord = lbrecord, count = count)
  }
  ##########################################################

  ################    Start simulation    ##################

  ##########################################################
  sigmoid = function(x) 1/(1+exp(-x))
  N_ = 300
  x = runif(N_, 0, 1)
  y = rbinom(length(x), 1, sigmoid(FUN(x)))

  fit = VB(y, x, rt0 = 0.01, st0 = 0.01, w0 = 1, J = J, tol = 1.0e-6)
  varphi = sqrt(2)*cos(outer(x,pi*(1:J)))
  unknown_f = varphi %*% fit$mutq

  xord = order(x)
  plot(x[xord], unknown_f[xord], col = 'purple', xlab = 'index', ylab = 'observed/fitted', main = 'Simulation result', ylim = c(-1, 1), type = 'l')
  curve(FUN, from = 0, to = 1, col = 'red', lty = 2, add = TRUE)
  legend("topright", legend = c('fitted values', 'true function'), col = c('purple', 'red'), lty = c(1, 2), bg = 'gray95')
  list(fit = fit, unknown_f = unknown_f, x = x)
}