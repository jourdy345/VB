source('cosineProbitMonotone.R')
simCosineProbitMonotone <- function(FUN, J, delta) {
  ##########################################################

  ################    Start simulation    ##################

  ##########################################################
  N_ <- 1000
  x <- runif(N_, 0, 1)
  W <- as.matrix(rep(1, times = N_))
  muBeta0 <- 0
  y <- rbinom(length(x),1,pnorm(FUN(x) + W + rnorm(N_)))
  SigmaBeta0 <- matrix(100,nrow=1,ncol=1)
  n <- length(x)
  muTheta0 <- rep(1, J+1)
  SigmaTheta0 <- diag(1, J+1)
  fit <- VB(y, x, W, J = J, delta = delta, r0s = 0.01, s0s = 0.01, r0t = 0.01, s0t = 0.01, w0 = 2, muBeta0 = muBeta0, SigmaBeta0 = SigmaBeta0, muTheta0, SigmaTheta0, sigma0=100, tol = 1.0e-06)
  
  varphi <- setVarPhi(x, J)
  unknownFunction <- rep(0, n)
  for (i in 1:n) {
    unknownFunction[i] <- sum(W[i] * fit$muBetaq) + delta * (sum(varphi[,,i] * tcrossprod(fit$muThetaq)) + sum(varphi[,,i] * fit$SigmaThetaq))
  }
  ord = order(x)
  plot(x[ord], unknownFunction[ord], col = 'purple', xlab = 'index', ylab = 'observed/fitted', ylim = c(0, 5), main = 'Simulation result', type = 'l')
  curve(FUN, from = 0, to = 1, col = 'red', lty = 2, add = TRUE)
  points(x,y,pch=19,col='grey')
  legend("topright", legend = c('fitted values', 'true function'), col = c('purple', 'red'), lty = c(1, 2), bg = 'gray95')
  list(fit = fit, y = y, x = x, W = W, unknownFunction = unknownFunction)
}


# myFun <- function(x) exp(0.2 * x) + 0.4*x - log(pi*x)
# myFun2 <- function(x) sech(1.127*pi*x)*log(9.46*x) + gamma(x) // monotonically decreasing

# ord <- order(VB_fit$x)
# plot(VB_fit$x[ord], VB_fit$unknownFunction[ord], col = 'purple', type = 'l')



# dOverp <- function(t) {
#   temp <- sqrt(2/pi)
#   temp/((temp * (-1/t + 1/t^3 - 1/t^5)) + 2 * exp(t^2/2))
# }