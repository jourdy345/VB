library(mvtnorm)    # for multivariate Gaussian distribution
library(matrixcalc) # for duplication matrix and elimination matrix
library(MCMCpack)   # for vech and xpnd
library(Matrix)     # for block-diagonal operator


# sample multivariate Gaussian with Cholesky factor of the precision matrix
rMNorm <- function(n,mu,C) {
  p <- length(mu)
  X <- matrix(rnorm(n*p),n,p)
  t(apply(X,1,function(x) solve(t(C),x)+mu))
}

logh <- function(n,p,y,X,beta,mu0,sig0,sigma2,A,B) {
  -0.5*n*log(2*pi*sigma2)-0.5/sigma2*sum((y-X%*%beta)*(y-X%*%beta))-0.5*determinant(sig0)$modulus[1]-0.5*p*log(2*pi)-
  0.5*sum((beta-mu0)*(solve(sig0,beta-mu0)))-(A+1)*log(sigma2)-B/sigma2+A*log(B)-lgamma(A)
}

logq <- function(p,beta,muq,C,sigma2,Aq,Bq) {
  -p*0.5*log(2*pi)+determinant(C)$modulus[1]-0.5*sum((crossprod(C,beta-muq))^2)-(Aq+1)*log(sigma2)-
  Bq/sigma2+Aq*log(Bq)-lgamma(Aq)
}


n <- 100
p <- 4

beta_true <- c(0.3, 10, 2, 6)
X <- matrix(rnorm(400), 100,4)
y <- as.vector(X%*%beta_true)+rnorm(100,sd=1.3)

sig0 <- diag(p)
mu0 <- rep(0,p)
A <- 1
B <- 1

# variational parameters
# sigq <- diag(p)
sigq <- crossprod(matrix(runif(p*p),p,p))
invsigq <- chol2inv(chol(sigq))
C <- t(chol(invsigq))
diag(C) <- log(diag(C))
muq <- rep(0,p)
Aq <- 1
Bq <- 1

L <- elimination.matrix(p)
D <- duplication.matrix(p)
Dmp <- ginv(D) # Moore-Penrose inverse


count <- 0
tau <- 1
kappa <- 0.8
lbold <- -Inf
lbnew <- 0
lb <- c()
S <- 5000

# for (i in 1:5000) {
#   count <- count+1
#   # rho <- 1/(count+tau)^kappa
#   rho <- 1/(count+1)^kappa
#   beta_sample <- rmvnorm(S,mean=muq,sigma=sigq)
#   sigma2_sample <- 1/rgamma(S,shape=Aq,rate=Bq)
#   # logh_v <- 0
#   # logq_v <- 0
#   muq_s <- rep(0,p)
#   sigq_s <- matrix(0,p,p)
#   Aq_s <- 0
#   Bq_s <- 0
#   for (j in 1:S) {
#     logh_v <- logh(n,p,y,X,beta_sample[j,],mu0,sig0,sigma2_sample[j],A,B)
#     logq_v <- logq(p,beta_sample[j,],muq,sigq,sigma2_sample[j],Aq,Bq)
#     f_v    <- logh_v-logq_v
#     muq_s   <- muq_s+invsigq%*%(beta_sample[j,]-muq)*f_v
#     sigq_s <- sigq_s-0.5*(invsigq+tcrossprod(beta_sample[j,]-muq))*f_v
#     Aq_s   <- Aq_s+(-log(sigma2_sample[j]+log(Bq)-digamma(Aq)))*f_v
#     Bq_s   <- Bq_s+(-1/sigma2_sample[j]+Aq/Bq)*f_v
#   }
#   muq <- muq+rho*muq_s/S
#   sigq <- sigq+rho*sigq_s/S
#   invsigq <- chol2inv(chol(sigq))
#   Aq <- Aq+rho*Aq_s/S
#   Bq <- Bq+rho*Bq_s/S
#   cat('Aq: ',Aq,', Bq: ',Bq,'\n')
#   # Lower bound
#   beta_sample <- rmvnorm(S,mean=muq,sigma=sigq)
#   sigma2_sample <- 1/rgamma(S,shape=Aq,rate=Bq)
#   lb_s <- 0
#   for (j in 1:S) {
#     logh_v <- logh(n,p,y,X,beta_sample[j,],mu0,sig0,sigma2_sample[j],A,B)
#     logq_v <- logq(p,beta_sample[j,],muq,sigq,sigma2_sample[j],Aq,Bq)
#     lb_s <- lb_s+logh_v-logq_v
#   }
#   lbnew <- lb_s/S
#   dif <- (lbnew-lbold)/lbold
#   lb <- c(lb,lbnew)
#   lbold <- lbnew
#   cat('lbnew: ',lbnew,', dif: ',dif,'\n')
# }
gauss_theta <- c(muq,vech(C))
invgam_theta <- c(Aq,Bq)
for (i in 1:5000) {
  count <- count+1
  rho <- 1/(count+1)^kappa
  beta_sample <- rMNorm(S,mu=muq,C=C)
  sigma2_sample <- 1/rgamma(S,shape=Aq,rate=Bq)
  gauss_grad <- 0
  invgam_grad <- 0
  for (j in 1:S) {
    logh_v <- logh(n,p,y,X,beta_sample[j,],mu0,sig0,sigma2_sample[j],A,B)
    logq_v <- logq(p,beta_sample[j,],muq,C,sigma2_sample[j],Aq,Bq)
    f_v    <- logh_v-logq_v
    gauss_grad <- gauss_grad+f_v*c(tcrossprod(C)%*%(beta_sample[j,]-muq),vech(1/diag(C)-tcrossprod(beta_sample[j,]-muq)%*%C))
    invgam_grad <- invgam_grad+f_v*c(log(Bq)-digamma(A)-log(sigma2_sample[j]),Aq/Bq-1/sigma2_sample[j])
  }
  gauss_grad <- gauss_grad/S
  invgam_grad <- invgam_grad/S
  gauss_invhess <- as.matrix(bdiag(solve(tcrossprod(C)),solve(2*L%*%(t(C)%x%diag(p))%*%D%*%Dmp%*%(sigq%x%sigq)%*%crossprod(Dmp,t(D))%*%tcrossprod(C%x%diag(p),L))))
  invgam_invhess <- matrix(c(Aq/(Bq*Bq),-Aq/Bq*(digamma(Aq)-digamma(Aq+1)),-Aq/Bq*(digamma(Aq)-digamma(Aq+1)),trigamma(Aq)),2,2)
  gauss_theta <- gauss_theta+rho*gauss_invhess%*%gauss_grad
  invgam_theta <- invgam_theta+rho*invgam_invhess%*%invgam_grad
  muq <- gauss_theta[1:p]
  C <- xpnd(gauss_theta[(p+1):length(gauss_theta)],p)*lower.tri(C,diag=TRUE)
  diag(C) <- exp(diag(C))
  print(muq)
  Aq <- invgam_theta[1]
  Bq <- invgam_theta[2]
  print(Aq)
  print(Bq)
  beta_sample <- rMNorm(S,mu=muq,C=C)
  sigma2_sample <- 1/rgamma(S,shape=Aq,rate=Bq)
  lb_s <- 0
  for (j in 1:S) {
    logh_v <- logh(n,p,y,X,beta_sample[j,],mu0,sig0,sigma2_sample[j],A,B)
    logq_v <- logq(p,beta_sample[j,],muq,C,sigma2_sample[j],Aq,Bq)
    lb_s <- lb_s+logh_v-logq_v
  }
  lbnew <- lb_s/S
  dif <- (lbnew-lbold)/abs(lbold)
  lb <- c(lb,lbnew)
  lbold <- lbnew
  cat('lbnew: ',lbnew,', dif: ',dif,'\n')
}