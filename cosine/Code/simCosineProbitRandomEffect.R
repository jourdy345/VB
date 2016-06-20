source('cosineProbitRandomEffect.R')

simCosineProbitRandomEffect = function(FUN, J, draw = TRUE, intercept = TRUE) {
  ##########################################################

  ################    Start simulation    ##################

  ##########################################################
  N_ = 1000
  x = runif(N_, 0, 1)
  W = rep(1, times = N_)
  muBeta0 = 0
  y = rbinom(length(x),1,pnorm(FUN(x) + W))
  SigmaBeta0 = matrix(1,nrow=1,ncol=1)

  # Need to be initialized
  muu0 = 
  Sigmau0 =
  
  fit = VB(y, x, W, muBeta0, SigmaBeta0, w0 = 1, r0s = 0.01, s0s = 0.01, r0t = 0.01, s0t = 0.01, muPsi0 = 1, muu0, Sigmau0, J = J, tol = 1.0e-06)
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

}