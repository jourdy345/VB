library(mvtnorm)    # for multivariate Gaussian distribution
library(matrixcalc) # for duplication matrix and elimination matrix
library(MCMCpack)   # for vech and xpnd
library(Matrix)     # for block-diagonal operator


# sample multivariate Gaussian with Cholesky factor of the precision matrix
rMNorm <- function(n,mu,tC) {
  p <- length(mu)
  X <- matrix(rnorm(n*p),n,p)
  t(apply(X,1,function(x) solve(tC,x)+mu))
}

logh <- function(n,p,y,X,beta,mu0,sig0,sigma2,A,B) {
  v <- -0.5*n*log(2*pi*sigma2)-0.5/sigma2*sum((y-X%*%beta)*(y-X%*%beta))-0.5*determinant(sig0)$modulus[1]-0.5*p*log(2*pi)-
  0.5*sum((beta-mu0)*(solve(sig0,beta-mu0)))-(A+1)*log(sigma2)-B/sigma2+A*log(B)-lgamma(A)
  v
}

logq <- function(p,beta,muq,C,sigma2,Aq,Bq) {
  v <- -p*0.5*log(2*pi)+determinant(C)$modulus[1]-0.5*sum((crossprod(C,beta-muq))^2)-(Aq+1)*log(sigma2)-
  Bq/sigma2+Aq*log(Bq)-lgamma(Aq)
  v
}

# XX: X^{t}X
# Xy: X^{t}y
# invsig0: Sigma_{0}^{-1}
# invsig0mu0: Sigma_{0}^{-1} * mu_{0}^{-1}
loghbeta_grad <- function(Xy,XX,beta,sigma2,invsig0,invsig0mu0) {
  as.vector((Xy-XX%*%beta)/sigma2+(invsig0mu0-invsig0%*%beta))
}

n <- 100
p <- 4

beta_true <- c(0.3, 10, 2, 6)
X <- matrix(rnorm(400), 100,4)
y <- as.vector(X%*%beta_true)+rnorm(100,sd=1.3)
XX <- crossprod(X)
Xy <- crossprod(X,y)

sig0 <- diag(p)
invsig0 <- solve(sig0)
mu0 <- rep(0,p)
invsig0mu0 <- invsig0%*%mu0
A <- 1
B <- 1

# variational parameters
# sigq <- diag(p)
sigq <- crossprod(matrix(rnorm(p*p),p,p))
invsigq <- chol2inv(chol(sigq))
tC <- chol(invsigq)
C <- t(tC)
Cold <- C
# diag(C) <- log(diag(C))
muq <- rep(0,p)
Aq <- 1
Bq <- 1

L <- elimination.matrix(p)
D <- duplication.matrix(p)
Dmp <- ginv(D) # Moore-Penrose inverse


count <- 0
tau <- 5
kappa <- 0.8
lbold <- -Inf
lb <- 0
lbc <- c()
S <- 25
K <- 5000
sigma2 <- 1.3*1.3
invgam_theta <- c(Aq,Bq)
for (i in 1:5000) {
  count <- count+1
  rho <- 1/(count+5)
  z <- rnorm(p)
  sigma2 <- 1/rgamma(1,shape=Aq,rate=Bq)
  # update Gaussian variational parameters
  beta <- solve(tC,z)+muq
  muq_grad <- loghbeta_grad(Xy,XX,beta,sigma2,invsig0,invsig0mu0)+as.vector(C%*%z)
  C_grad <- -tcrossprod(solve(tC,z),solve(C,muq_grad))
  diag(C_grad) <- diag(C_grad)*diag(C)
  muq <- muq+rho*muq_grad
  C <- C+rho*C_grad
  # diag(C) <- exp(diag(C))
  # C <- C*lower.tri(C,diag=TRUE)
  tC <- t(C)
  sigq <- tcrossprod(C)
  if (any(eigen(sigq)$values<0)) {
    C <- Cold
  } else {
    Cold <- C
  }
  # update inverse-gamma variational parameters
  sigma2_sample <- 1/rgamma(S,shape=Aq,rate=Bq)
  cntrl_var <- sapply(sigma2_sample,function(x) c(log(Bq)-digamma(Aq)-log(x),Aq/Bq-1/x))
  m <- rowMeans(cntrl_var)
  h_v <- sapply(sigma2_sample,function(x) logq(p,beta,muq,C,x,Aq,Bq)-logh(n,p,y,X,beta,mu0,sig0,x,A,B))
  invgam_grad <- colSums(apply(apply(cntrl_var,2,function(x) x-m),1,function(x) x*h_v))/S
  invgam_hess <- tcrossprod(apply(cntrl_var,2,function(x) x-m))/(S-1)
  invgam_theta <- invgam_theta+rho*solve(invgam_hess,invgam_grad)
  print(invgam_theta)
  Aq <- invgam_theta[1]
  Bq <- invgam_theta[2]
  # cat('Aq: ',Aq,', Bq: ',Bq,'\n')
  beta_sample <- rMNorm(K,mu=muq,tC=tC)
  sigma2_sample <- 1/rgamma(K,shape=Aq,rate=Bq)
  lb <- 0
  for (k in 1:K) {
    lb <- lb+(logh(n,p,y,X,beta_sample[k,],mu0,sig0,sigma2_sample[k],A,B)-logq(p,beta_sample[k,],muq,C,sigma2_sample[k],Aq,Bq))
  }
  lb <- lb/K
  dif <- (lb-lbold)/abs(lbold)
  lbold <- lb
  cat('lb: ',lb,', dif: ',dif,'\n')
}

