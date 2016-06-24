source('cosineProbitMonotone.R')
simCosineProbitMonotone <- function(FUN, J) {
  ##########################################################

  ################    Start simulation    ##################

  ##########################################################
  N_ <- 1000
  x <- runif(N_, 0, 1)
  W <- as.matrix(rep(1, times = N_))
  muBeta0 <- 0
  y <- rbinom(length(x),1,pnorm(FUN(x) + W + rnorm(N_)))
  SigmaBeta0 <- matrix(1,nrow=1,ncol=1)
  n <- length(x)
  delta <- -1
  muTheta0 <- rep(0, J+1)
  SigmaTheta0 <- diag(1, J+1)
  fit <- VB(y, x, W, J = J, delta = delta, r0s = 0.01, s0s = 0.01, r0t = 0.01, s0t = 0.01, w0 = 1, muBeta0 = muBeta0, SigmaBeta0 = SigmaBeta0, muTheta0, SigmaTheta0, sigma0=1, tol = 1.0e-04)
  
  varphi <- setVarPhi(x, J)
  unknownFunction <- rep(0, n)
  for (i in 1:n) {
    unknownFunction[i] <- sum(W[i] * fit$muBetaq) + delta * sum(varphi[,,i] * tcrossprod(fit$muThetaq))
  }
  ord = order(x)
  plot(x[ord], unknownFunction[ord], col = 'purple', xlab = 'index', ylab = 'observed/fitted', main = 'Simulation result', ylim = c(-5, 5), type = 'l')
  curve(FUN, from = 0, to = 1, col = 'red', lty = 2, add = TRUE)
  legend("topright", legend = c('fitted values', 'true function'), col = c('purple', 'red'), lty = c(1, 2), bg = 'gray95')
  list(fit = fit, y = y, x = x)
}


# myFun <- function(x) exp(0.2 * x) + 0.4*x - log(pi*x)