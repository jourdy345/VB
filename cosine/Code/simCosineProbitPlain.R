source('cosineProbitPlain.R')
simCosineProbitPlain <- function(FUN, J = 20, tol = 1.0e-06) {
  N_ <- 1000
  x <- runif(N_, 0, 1)
  W <- rep(1, times = N_)
  muBeta0 <- 0
  y <- rbinom(length(x),1,pnorm(FUN(x) + W))
  SigmaBeta0 <- matrix(1,nrow=1,ncol=1)
  fit <- VB(x, y, W, muBeta0 = muBeta0, SigmaBeta0 = SigmaBeta0, w0 = 1, r0s = 0.01, s0s = 0.01, r0t = 0.01, s0t = 0.01, J = J, tol = 1.0e-06)
  
  varphi <- sqrt(2) * cos(outer(x, pi * (1:J)))
  unknownFunction <- W %*% fit$muBetaq + varphi %*% fit$muThetaq
  ord <- order(x)
  plot(x[ord], unknownFunction[ord], col = 'purple', xlab = 'index', ylab = 'observed/fitted', main = 'Simulation result', ylim = c(-10, 10), type = 'l')
  curve(FUN, from = 0, to = 1, col = 'red', lty = 2, add = TRUE)
  legend('topright', legend = c('fitted values', 'true function'), col = c('purple', 'red'), lty = c(1, 2), bg = 'gray95')
  list(fit = fit, y = y, x = x, varphi = varphi)
}

# myFun <- function(x) sin(2 * pi * x) * log(x / 5)