library(mvtnorm)
library(MCMCpack)

log_kernel_alpha <- function(y,b,lam,alpha,sig.a2) {
  # v <- prod(exp(-0.5*(lam+exp(alpha)*b))*exp(-0.5*exp(-lam-exp(alpha)*b)*y*y)*exp(-0.5/sig.a2*alpha*alpha))
  v <- sum(-0.5*(exp(alpha)*b)-0.5*exp(-lam-exp(alpha)*b)*y*y-0.5/sig.a2*alpha*alpha)
  v
}

log_kernel_b1 <- function(y,b,lam,alpha,phi) {
  # v <- exp(-0.5*(lam+exp(alpha)*b[1]))*exp(-0.5*exp(-lam-exp(alpha)*b[1])*y[1]^2)*exp(-0.5*(1-phi^2)*b[1]^2)*exp(-0.5*(b[2]-phi*b[1])^2)
  v <- -0.5*(exp(alpha)*b[1])-0.5*exp(-lam-exp(alpha)*b[1])*y[1]^2-0.5*(1-phi^2)*b[1]^2-0.5*(b[2]-phi*b[1])^2
  v
}

log_kernel_bt <- function(y,b,lam,alpha,phi,n) {
  # v <- exp(-0.5*(lam+exp(alpha)*b[2:(n-1)]))*exp(-0.5*exp(-lam-exp(alpha)*b[2:(n-1)])*y[2:(n-1)]^2)*exp(-0.5*(b[2:(n-1)]-phi*b[1:(n-2)])^2-0.5*(b[3:n]-phi*b[2:(n-1)]^2))
  v <- -0.5*(exp(alpha)*b[2:(n-1)])-0.5*exp(-lam-exp(alpha)*b[2:(n-1)])*y[2:(n-1)]^2-0.5*(b[2:(n-1)]-phi*b[1:(n-2)])^2-0.5*(b[3:n]-phi*b[2:(n-1)]^2)
  v
}

log_kernel_bn <- function(y,b,lam,alpha,phi,n) {
  # v <- exp(-0.5*(lam+exp(alpha)*b[n]))*exp(-0.5*exp(-lam-exp(alpha)*b[n])*y[n]^2)*exp(-0.5*(b[n]-phi*b[n-1])^2)
  v <- -0.5*(exp(alpha)*b[n])-0.5*exp(-lam-exp(alpha)*b[n])*y[n]^2-0.5*(b[n]-phi*b[n-1])^2
  v
}

log_kernel_lam <- function(y,b,lam,alpha,sig.l2) {
  # v <- prod(exp(-0.5*(lam+exp(alpha)*b))*exp(-0.5*exp(-lam-exp(alpha)*b)*y*y))*exp(-0.5/sig.l2*lam*lam)
  v <- -0.5*n*lam-0.5*sum(exp(-lam-exp(alpha)*b)*y*y)-0.5/sig.l2*lam*lam
  v
}

log_kernel_psi <- function(b,psi,n) {
  # v <- sqrt(1-(exp(psi)/(1+exp(psi)))^2)*exp(-0.5*(1-(exp(psi)/(1+exp(psi)))^2)*b[1]^2)*prod(exp(-0.5*(b[2:n]-(exp(psi)/(1+exp(psi)))*b[1:(n-1)])^2))*exp(-0.5/sig.p2*psi*psi)
  v <- 0.5*log(1-(exp(psi)/(1+exp(psi)))^2)-0.5*(1-(exp(psi)/(1+exp(psi)))^2)*b[1]^2-0.5*sum((b[2:n]-(exp(psi)/(1+exp(psi)))*b[1:(n-1)])^2)-0.5/sig.p2*psi*psi
  v
}


n <- 60
d <- n+3
psi_true <- 0.3
phi_true <- 1/(1+exp(-psi_true))
lam_true <- 1.24
sig_true <- 2.73
b1_true <- rnorm(1,0,sqrt(1/(1-phi_true*phi_true)))
b_true <- rep(0,n)
b_true[1] <- b1_true
y <- rep(0,n)
y[1] <- rnorm(1,0,exp(0.5*(lam_true+sig_true*b_true[1])))
for (t in 2:n) {
  b_true[t] <- rnorm(1,phi_true*b_true[t-1],1)
  y[t] <- rnorm(1,0,exp(0.5*(lam_true+sig_true*b_true[t])))
}

sig.a2 <- sig.l2 <- sig.p2 <- 5
alpha_true <- log(sig_true)


burnin <- 100000
thinin <- 10
nSample <- 30000

alpha <- 1
b <- rnorm(n)
lam <- psi <- 1
sig <- exp(lam)
phi <- exp(psi)/(1+exp(psi))

for (i in 1:burnin) {
  alpha_prop <- rnorm(1)
  b_prop     <- rnorm(n)
  lam_prop <- rnorm(1)
  psi_prop <- rnorm(1)
  sig_prop <- exp(lam_prop)
  phi_prop <- exp(psi_prop)/(1+exp(psi_prop))
  accept_alpha <- exp(log_kernel_alpha(y,b,lam,alpha_prop,sig.a2)-log_kernel_alpha(y,b,lam,alpha,sig.a2)+dnorm(alpha,log=TRUE)-dnorm(alpha_prop,log=TRUE))
  accept_b1    <- exp(log_kernel_b1(y,b_prop,lam,alpha,phi)-log_kernel_b1(y,b,lam,alpha,phi)+dnorm(b[1],log=TRUE)-dnorm(b_prop[1],log=TRUE))
  accept_bt    <- exp(log_kernel_bt(y,b_prop,lam,alpha,phi,n)-log_kernel_bt(y,b,lam,alpha,phi,n)+dnorm(b[2:(n-1)],log=TRUE)-dnorm(b_prop[2:(n-1)],log=TRUE))
  accept_bn    <- exp(log_kernel_bn(y,b_prop,lam,alpha,phi,n)-log_kernel_bn(y,b,lam,alpha,phi,n)+dnorm(b[n],log=TRUE)-dnorm(b_prop[n],log=TRUE))
  accept_lam   <- exp(log_kernel_lam(y,b,lam_prop,alpha,sig.l2)-log_kernel_lam(y,b,lam,alpha,sig.l2)+dnorm(lam,log=TRUE)-dnorm(lam_prop,log=TRUE))
  accept_psi   <- exp(log_kernel_psi(b,psi_prop,n)-log_kernel_psi(b,psi,n)+dnorm(psi_prop,log=TRUE)-dnorm(psi,log=TRUE))

  accept_b     <- c(accept_b1,accept_bt,accept_bn)
  binary_b     <- 1*(runif(n) < accept_b)
  if (runif(1)<accept_alpha) alpha <- alpha_prop
  b <- b_prop*binary_b+b*(1-binary_b)
  if (runif(1)<accept_lam) lam <- lam_prop
  if (runif(1)<accept_psi) psi <- psi_prop
}

theta_storage <- matrix(0,n+3,nSample)
for (i in 1:nSample) {
  for (j in 1:thinin) {
    alpha_prop <- rnorm(1)
    b_prop     <- rnorm(n)
    lam_prop <- rnorm(1)
    psi_prop <- rnorm(1)
    sig_prop <- exp(lam_prop)
    phi_prop <- exp(psi_prop)/(1+exp(psi_prop))
    accept_alpha <- exp(log_kernel_alpha(y,b,lam,alpha_prop,sig.a2)-log_kernel_alpha(y,b,lam,alpha,sig.a2)+dnorm(alpha,log=TRUE)-dnorm(alpha_prop,log=TRUE))
    accept_b1    <- exp(log_kernel_b1(y,b_prop,lam,alpha,phi)-log_kernel_b1(y,b,lam,alpha,phi)+dnorm(b[1],log=TRUE)-dnorm(b_prop[1],log=TRUE))
    accept_bt    <- exp(log_kernel_bt(y,b_prop,lam,alpha,phi,n)-log_kernel_bt(y,b,lam,alpha,phi,n)+dnorm(b[2:(n-1)],log=TRUE)-dnorm(b_prop[2:(n-1)],log=TRUE))
    accept_bn    <- exp(log_kernel_bn(y,b_prop,lam,alpha,phi,n)-log_kernel_bn(y,b,lam,alpha,phi,n)+dnorm(b[n],log=TRUE)-dnorm(b_prop[n],log=TRUE))
    accept_lam   <- exp(log_kernel_lam(y,b,lam_prop,alpha,sig.l2)-log_kernel_lam(y,b,lam,alpha,sig.l2)+dnorm(lam,log=TRUE)-dnorm(lam_prop,log=TRUE))
    accept_psi   <- exp(log_kernel_psi(b,psi_prop,n)-log_kernel_psi(b,psi,n)+dnorm(psi_prop,log=TRUE)-dnorm(psi,log=TRUE))

    accept_b     <- c(accept_b1,accept_bt,accept_bn)
    binary_b     <- 1*(runif(n) < accept_b)
    if (runif(1)<accept_alpha) alpha <- alpha_prop
    b <- b_prop*binary_b+b*(1-binary_b)
    if (runif(1)<accept_lam) lam <- lam_prop
    if (runif(1)<accept_psi) psi <- psi_prop
  }
  theta_storage[,i] <- c(b,alpha,lam,psi)
}

MCMCmeans <- rowMeans(theta_storage)
MCMCest   <- cbind(c(b_true,alpha_true,lam_true,psi_true),MCMCmeans)