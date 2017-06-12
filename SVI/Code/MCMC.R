library(mvtnorm)
library(msm)        # for half-normal random variate generation
library(matrixcalc) # for duplication matrix and elimination matrix
library(MCMCpack)   # for vech and xpnd
library(Matrix)     # for block-diagonal operator


n <- 100
p <- 4

X <- matrix(rnorm(n*p),n,p)*10
beta_true <- c(1.3,2.5,10,-0.24)
y <- X%*%beta_true+rnorm(n)*1.3


XX <- crossprod(X)
Xy <- crossprod(X,y)
sig <- diag(p)
invsig <- chol2inv(chol(sig))
mu <- rep(0,p)
invsig_mu <- invsig%*%mu
a <- 2

joint <- function(y,X,beta,sigma2,a,n,p) {
  exp(-0.5/sigma2*sum((y-X%*%beta)^2)-sigma2^2/(2*a^2))/sqrt(sigma2)^n
}

beta <- rnorm(p)
sigma2 <- rexp(1)

# infer with MCMC
for (i in 1:10000) {
  sigq <- solve(1/sigma2*XX+invsig)
  muq  <- sigq%*%(1/sigma2*Xy+invsig_mu)
  beta <- c(rmvnorm(1,muq,sigq))
  sig.prop <- rexp(1)
  a1 <- joint(y,X,beta,sig.prop,a,n,p)/joint(y,X,beta,sigma2,a,n,p)
  a2 <- dexp(sigma2)/dexp(sig.prop)
  if (runif(1)<(a1*a2)) sigma2 <- sig.prop
}

thetaContainer <- matrix(0,p+1,10000)
for (i in 1:10000) {
  for (j in 1:4) {
    sigq <- solve(1/sigma2*XX+invsig)
    muq  <- sigq%*%(1/sigma2*Xy+invsig_mu)
    beta <- c(rmvnorm(1,muq,sigq))
    sig.prop <- rexp(1)
    a1 <- joint(y,X,beta,sig.prop,a,n,p)/joint(y,X,beta,sigma2,a,n,p)
    a2 <- dexp(sigma2)/dexp(sig.prop)
    if (runif(1)<(a1*a2)) sigma2 <- sig.prop
  }
  thetaContainer[,i] <- c(beta,sigma2)
  cat('i: ',i,'\n')
}

# variational inference

# logh <- function(y,X,theta,mu,invsig,a,n,p) {
#   beta <- theta[1:p]
#   sigma2 <- theta[p+1]
#   # v <- dmvnorm(x=y,mean=c(X%*%beta),sigma=diag(sigma2,n),log=TRUE)+dmvnorm(x=beta,mean=mu,sigma=sig,log=TRUE)+log(sqrt(2/pi)/a*exp(-sigma2^2/(2*a^2)))
#   -0.5/sigma2*sum((y-X%*%beta)^2)-0.5*sum((beta-mu)*(invsig%*%(beta-mu)))-sigma2^2/(2*a^2)-n/2*log(sigma2)
# }

# dtheta <- function(XX,Xy,y,X,theta,mu,invsig,a,n,p) {
#   beta <- theta[1:p]
#   sigma2 <- theta[p+1]
#   dbeta <- c(1/sigma2*(Xy-XX%*%beta)-invsig%*%(beta-mu))
#   dsigma2 <- 1/(2*sigma2^2)*sum((y-X%*%beta)^2)-sigma2/a^2-n/(2*sigma2)
#   c(dbeta,dsigma2)
# }

# muq <- rep(0,p+1)
# C   <- diag(p+1)
# nIter <- 5000
# lbc <- c()
# for (i in 1:nIter) {
#   rho <- 1/(5+i)
#   z <- rnorm(p+1)
#   theta <- c(C%*%z)+muq
#   dt <- dtheta(XX,Xy,y,X,theta,mu,invsig,a,n,p)
#   muq <- muq+rho*dt
#   C <- C+rho*(tcrossprod(dt,z)+diag(1/diag(C)))
#   lb <- logh(y,X,theta,mu,invsig,a,n,p)+determinant(C)$modulus[1]
#   lbc <- c(lbc,lb)
#   cat('lb: ',lb,'\n')
# }

logh <- function(y,X,beta,mu0,sig0,tau,A,n) {
  -n/2*tau-exp(-tau)/2*sum((y-X%*%beta)^2)+tau-0.5*sum((beta-mu0)*solve(sig0,beta-mu0))-exp(2*tau)/(2*A^2)
}

logq <- function(beta,muq,sigq,tau,mu.t,sig.t2) {
  -0.5*determinant(sigq)$modulus[1]-0.5*sum((beta-muq)*solve(sigq,beta-muq))-0.5*log(sig.t2)-0.5/sig.t2*(tau-mu.t)^2
}

grad_muq_fun <- function(beta,muq,sigq) {
  -solve(sigq,(muq-beta))
}
grad_sigq_fun <- function(beta,muq,sigq) {
  -0.5*solve(sigq)-0.5*tcrossprod(beta-muq)
}

grad_mut_fun <- function(tau,mu.t,sig.t2) {
  (tau-mu.t)/sig.t2
}

grad_sig.t2_fun <- function(tau,mu.t,sig.t2) {
  -0.5/sig.t2+0.5/sig.t2^2*(tau-mu.t)^2
}





n <- 100
p <- 4

X <- matrix(rnorm(n*p),n,p)*10
beta_true <- c(1.3,2.5,10,-0.24)
y <- X%*%beta_true+rnorm(n)*1.3


XX <- crossprod(X)
Xy <- crossprod(X,y)
sig0 <- diag(p)
invsig <- chol2inv(chol(sig))
mu0 <- rep(0,p)
invsig_mu <- invsig%*%mu
A <- 2

# variational parameters
muq <- rep(0,p)
sigq <- diag(p)
mu.t <- 0
sig.t2 <- 2


S <- 40
nIter <- 5000
LB <- 0
LBC <- c()

for (i in 1:nIter) {
  muq_old <- muq
  sigq_old <- sigq
  mu.t.old <- mu.t
  sig.t2.old <- sig.t2
  rho <- 1/(5+i)

  # generate MC samples
  beta_sample <- rmvnorm(S,muq,sigq)
  tau_sample  <- rnorm(S,mu.t,sig.t2)

  # compute control variates
  muq.m <- rowMeans(apply(beta_sample,1,function(beta) grad_muq_fun(beta,muq,sigq)))
  sigq.m <- matrix(0,p,p)
  mu.t.m <- 0
  sig.t2.m <- 0
  for (j in 1:S) {
    beta <- beta_sample[j,]
    tau  <- tau_sample[j]
    sigq.m <- sigq.m+grad_muq_fun(beta,muq,sigq)
    mu.t.m <- mu.t.m+grad_mut_fun(tau,mu.t,sig.t2)
    sig.t2.m <- sig.t2.m+grad_sig.t2_fun(tau,mu.t,sig.t2)
  }
  sigq.m <- sigq.m/S
  mu.t.m <- mu.t.m/S
  sig.t2.m <- sig.t2.m/S

  # compute gradients
  grad_muq <- 0
  grad_siq <- 0
  grad_mu.t <- 0
  grad_sig.t2 <- 0

  for (j in 1:S) {
    beta <- beta_sample[j,]
    tau  <- tau_sample[j]
    logh_v <- logh(y,X,beta,mu0,sig0,tau,A,n)
    logq_v <- logq(beta,muq,sigq,tau,mu.t,sig.t2)
    f_v    <- logh_v-logq_v
    LB     <- LB+f_v
    grad_muq <- grad_muq+(grad_muq_fun(beta,muq,sigq)-muq.m)*f_v
    grad_sigq <- grad_sigq+(grad_sigq_fun(beta,muq,sigq)-sigq.m)*f_v
    grad_mu.t <- grad_mu.t+(grad_mut_fun(tau,mu.t,sig.t2)-mu.t.m)*f_v
    grad_sig.t2 <- grad_sig.t2+(grad_sig.t2_fun(tau,mu.t,sig.t2)-sig.t2.m)*f_v
  }
  grad_muq <- grad_muq/(S-1)
  grad_sigq <- grad_sigq/(S-1)
  grad_mu.t <- grad_mu.t/(S-1)
  grad_sig.t2 <- grad_sig.t2/(S-1)
  LB <- LB/S

  # update variational parameters
  muq <- (1-rho)*muq+rho*grad_muq
  sigq <- (1-rho)*sigq+rho*grad_sigq
  mu.t <- (1-rho)*mu.t+rho*grad_mu.t
  sig.t2 <- (1-rho)*sig.t2+rho*grad_sig.t2
  if (any(eigen(sigq)$values < 0)) {
    sigq <- sigq_old
    muq <- muq_old
  }
  if (sig.t2 < 0) {
    sig.t2 <- sig.t2.old
    mu.t   <- mu.t.old
  }
  LBC <- c(LBC,LBC)
  cat('LB=',LB,'\n')
}











