source('../BSAR/Gbsar_v2.R')
source('../BSAR/Gbsar_probit_v2.R')
source('../BSAR/Gbsar_fun_Module_v2.R')
source('cosineProbitPlain.R')
library(BSAM)
library(gam)
simCosineProbitPlain <- function(FUN, J = 20, tol = 1.0e-06) {
  N_ <- 1000
  x <- runif(N_, 0, 1)
  W <- rep(1, times = N_)
  muBeta0 <- 0
  y <- rbinom(length(x),1,pnorm(FUN(x) + W))
  SigmaBeta0 <- matrix(1,nrow=1,ncol=1)
  fit <- VB(x, y, W, muBeta0 = muBeta0, SigmaBeta0 = SigmaBeta0, w0 = 1, r0t = 0.001, s0t = 0.001, J = J, tol = 1.0e-06)
  varphi <- sqrt(2) * cos(outer(x, pi * (1:J)))
  unknownFunction <- W %*% fit$muBetaq + varphi %*% fit$muThetaq
  # ord <- order(x)
  # plot(x[ord], unknownFunction[ord], col = 'purple', xlab = 'index', ylab = 'observed/fitted', main = 'Simulation result', ylim = c(-10, 10), type = 'l')
  # curve(FUN, from = 0, to = 1, col = 'red', lty = 2, add = TRUE)
  # legend('topright', legend = c('fitted values', 'true function'), col = c('purple', 'red'), lty = c(1, 2), bg = 'gray95')
  list(fit = fit, y = y, x = x, varphi = varphi, W = W, unknownFunction = unknownFunction)
}

### Comparison with prof. Seong-il Jo's
# par(mfrow = c(1, 2))
# myFun <- function(x) sin(2 * pi * x) * log(x / 5)
# myFun <- function(x) 2 * x + 1

# fit_VB <- simCosineProbitPlain(myFun, 20)
# fout_MCMC <- gbsar(y = fit_VB$y, w = as.matrix(fit_VB$W), x = fit_VB$x, family = 'bernoulli', link = 'probit', nbasis = 20, shape = 'Free')

# fit_MCMC <- fitted(fout_MCMC)
# o <- order(fit_VB$x)
# ylim <- range(c(myFun(fit_VB$x[o]), fit_MCMC$fxobs$mean, fit_VB$unknownFunction), na.rm = TRUE)
# # plot(fit_VB$x[o], fit_VB$unknownFunction[o], type = 'l', ylim = ylim, lwd = 2, ylab = '', xlab = '', main = 'sin(2*pi*x)*log(x/5)')
# plot(fit_VB$x[o], fit_VB$unknownFunction[o], type = 'l', ylim = ylim, lwd = 2, ylab = '', xlab = '', main = '2*x+1')
# lines(fit_VB$x[o], fit_MCMC$fxobs$mean[o], lwd = 2, col = 2, lty = 2)
# curve(myFun, from = 0, to = 1, col = 3, lwd = 2, lty = 3, add = TRUE)
# # points(fit_VB$x, fit_VB$y, pch=19, col='grey')
# rug(fit_VB$x)
# legend('topright', legend = c('fit_VB', 'fit_MCMC', 'True curve'), col = c(1, 2, 3), lty = c(1, 2, 3), bg = 'gray95')


### Comparison with prof. Seong-il Jo's
# myFun2 <- function(x) pi * x - 3 * sin(pi * x)
# myFun2 <- function(x) 1.5*sin(pi*x)
# myFun2 <- function(x) exp(2*x^2)-3
# fit_VB2 <- simCosineProbitPlain(myFun2, 30)
# fout_MCMC2 <- gbsar(y = fit_VB2$y, w = as.matrix(fit_VB2$W), x = fit_VB2$x, family = 'bernoulli', link = 'probit', nbasis = 30, shape = 'Free')
# fit_MCMC2 <- fitted(fout_MCMC2)
# o2 <- order(fit_VB2$x)
# ylim2 <- range(c(myFun2(fit_VB2$x[o2]), fit_MCMC2$fxobs$mean, fit_VB2$unknownFunction), na.rm = TRUE)
# plot(fit_VB2$x[o2], fit_VB2$unknownFunction[o2], type = 'l', ylim = ylim2, lwd = 2, ylab = '', xlab = '', main = '1.5*sin(pi*x)')
# lines(fit_VB2$x[o2], fit_MCMC2$fxobs$mean[o2], lwd = 2, col = 2, lty = 2)
# curve(myFun2, from = 0, to = 1, col = 3, lwd = 2, lty = 3, add = TRUE)
# rug(fit_VB2$x)
# legend('topright', legend = c('fit_VB', 'fit_MCMC', 'True curve'), col = c(1, 2, 3), lty = c(1, 2, 3), bg = 'gray95')

### Comparison with Dason Lee's
# myFun3 <- function(x) exp(2*x^2)-3
# myFun3 <- function(x) 2*x -3
myFun3 <- function(x) sin(2 * pi * x) * log(x / 5)
# myFun3 <- function(x) exp(2*x^2) - 3
nblow0=1000 # B0 iteration adaptive MCMC for shape-restricted model
nblow=1000 # Burn in
smcmc=1000 # Save
nskip=10 # Thin
ndisp=1000 # Display number of iteration
mcmc=list(nblow0=nblow0,nblow=nblow,smcmc=smcmc,nskip=nskip,ndisp=ndisp)

set.seed(123)
fit_VB3 <- simCosineProbitPlain(myFun3, 20)
GAM_dataset3 <- data.frame(y = fit_VB3$y, x = fit_VB3$x)
fit_Gbsar <- Gbsar(fit_VB3$y, xobs = fit_VB3$x, family = 'bernoulli', link = 'probit', method = 'DataAug', fmodel = 1, fpm = 1, nbasis = 20, mcmc = mcmc)
fit_GAM <- gam(y ~ s(x), family = binomial(link = probit), GAM_dataset3)
o3 <- order(fit_VB3$x)
ylim3 <- range(c(myFun3(fit_VB3$x[o3]), fit_VB3$unknownFunction, fit_Gbsar$post.est$muhatm[o3]), na.rm = TRUE)
plot(fit_VB3$x[o3], fit_VB3$unknownFunction[o3], type = 'l', ylim = ylim3, lwd = 2, ylab = '', xlab = '', main = 'sin(2 * pi * x) * log(x / 5)')
lines(fit_VB3$x[o3], fit_Gbsar$post.est$muhatm[o3], lwd = 2, col = 2, lty = 2)
lines(fit_VB3$x[o3],qnorm(fitted(fit_GAM)[o3]),col='darkgreen',lwd=2,lty=6)
lines(fit_VB3$x[o3], myFun3(fit_VB3$x[o3]) + 1, col = 3, lwd = 2, lty = 3)
# curve(myFun3, from = 0, to = 1, col = 3, lwd = 2, lty = 3, add = TRUE)
rug(fit_VB3$x)
legend('topright', legend = c('VB', 'MCMC', 'GAM', 'True curve'), col = c(1, 2, 'darkgreen', 3), lty = c(1, 2, 6, 3), bg = 'gray95')

