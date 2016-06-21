source('cosineProbitRandomEffect.R')

simCosineProbitRandomEffect = function(FUN, J, draw = TRUE, intercept = TRUE) {
  ##########################################################

  ################    Start simulation    ##################

  ##########################################################
  N_ = 1000
  x = runif(N_, 0, 1)
  W = as.matrix(rep(1, times = N_))
  Z = as.matrix(rep(1, times = N_))
  muBeta0 = 0
  y = rbinom(length(x),1,pnorm(FUN(x) + W + Z))
  SigmaBeta0 = matrix(1,nrow=1,ncol=1)

  # Need to be initialized
  muu0 = 0
  Sigmau0 = matrix(1, nrow = 1, ncol = 1)
  
  fit = VB(y, x, W, Z, muBeta0 = muBeta0, SigmaBeta0 = SigmaBeta0, w0 = 1, r0s = 0.01, s0s = 0.01, r0t = 0.01, s0t = 0.01, muPsi0 = 1, muu0 = muu0, Sigmau0 = Sigmau0, J = J, tol = 1.0e-04)
  categorize = function(x) {
    y = 0
    if (x < 0) {
      y = 0
    } else {
      y = 1
    }
    y
  }

  varphi = sqrt(2)*cos(outer(x,pi*(1:J)))
  unknownFunction = varphi %*% fit$muThetaq + W %*% fit$SigmaBetaq + Z %*% fit$muuq
  ord = order(x)
  plot(x[ord], unknownFunction[ord], col = 'purple', xlab = 'index', ylab = 'observed/fitted', main = 'Simulation result', ylim = c(-5, 5), type = 'l')
  curve(FUN, from = 0, to = 1, col = 'red', lty = 2, add = TRUE)
  legend("topright", legend = c('fitted values', 'true function'), col = c('purple', 'red'), lty = c(1, 2), bg = 'gray95')
  list(fit = fit, y = y, x = x)
}


