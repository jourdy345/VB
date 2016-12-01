###################################################
# SVI for linear model with control variate       #
# Salimans, T., & Knowles, D. A. (2014).          #
# On using control variates with stochastic       #
# approximation for variational Bayes and         #
# its connection to stochastic linear regression. #
# arXiv preprint arXiv:1401.1022.                 #
###################################################

library(mvtnorm)    # for multivariate Gaussian distribution
library(matrixcalc) # for duplication matrix and elimination matrix
library(MCMCpack)   # for vech and xpnd
library(Matrix)     # for block-diagonal operator



# q(beta)   ~ N(muq,sigq)
# q(sigma2) ~ InvGam(Aq,Bq)



# log-joint density (log(likelihood * prior))
logh <- function(n,p,y,X,beta,mu0,sig0,sigma2,A,B) {
  -0.5*n*log(2*pi*sigma2)-0.5/sigma2*sum((y-X%*%beta)*(y-X%*%beta))-0.5*determinant(sig0)$modulus[1]-0.5*p*log(2*pi)-
  0.5*sum((beta-mu0)*(solve(sig0,beta-mu0)))-(A+1)*log(sigma2)-B/sigma2+A*log(B)-lgamma(A)
}

# log-variational posterior
logq <- function(p,beta,muq,invsigq,sigma2,Aq,Bq) {
  -p*0.5*log(2*pi)+0.5*determinant(invsigq)$modulus[1]-0.5*sum((beta-muq)*(invsigq%*%(beta-muq)))-(Aq+1)*log(sigma2)-
  Bq/sigma2+Aq*log(Bq)-lgamma(Aq)
}

# generate simulation data
n <- 100
p <- 4
beta_true <- c(0.3, 10, 2, 6)
X <- matrix(rnorm(400), 100,4)
y <- as.vector(X%*%beta_true)+rnorm(100,sd=1.3)

# initialize hyperparameters for prior distributions
sig0 <- diag(p)
mu0 <- rep(0,p)
A <- 1
B <- 1

# variational parameters
sigq <- crossprod(matrix(runif(p*p),p,p))
invsigq <- chol2inv(chol(sigq))
muq <- rep(0,p)
Aq <- 1
Bq <- 1


L <- elimination.matrix(p)
D <- duplication.matrix(p)
D_plus <- ginv(D) # Moore-Penrose inverse
count <- 0
tau <- 1
kappa <- 0.8
lbold <- -Inf
lb <- 0
lbc <- c()
S <- 40 # number of Monte-Carlo samples
T <- 4000 # number of iterations

# natural parameters
gauss_lambda1 <- c(solve(sigq,muq))
gauss_lambda2 <- -0.5*crossprod(D,as.vector(invsigq))
gauss_lambda <- c(gauss_lambda1,gauss_lambda2)
invgam_lambda <- c(Aq,Bq)

muq_storage <- matrix(0,4,T)
sigma2_storage <- matrix(0,2,T)

# start inference
for (i in 1:T) {
  gauss_lambda_old  <- gauss_lambda
  muq_old           <- muq
  sigq_old          <- sigq
  invgam_lambda_old <- invgam_lambda
  Aq_old            <- Aq
  Bq_old            <- Bq
  rho   <- 1/(5+i)
  M_mat <- 2*D_plus%*%(muq%x%diag(p))
  S_mat <- 2*D_plus%*%tcrossprod((sigq%x%sigq),D_plus)
  aux   <- solve(S_mat,M_mat)
  gauss_invfisherinformation <- rbind(cbind(invsigq+crossprod(M_mat,aux),-t(aux)),cbind(-aux,solve(S_mat)))
  invgam_fisherinformation   <- matrix(c(trigamma(Aq),-1/Bq,-1/Bq,Aq/(Bq*Bq)),2,2)


  beta_sample   <- rmvnorm(S,mean=muq,sigma=sigq)
  sigma2_sample <- 1/rgamma(S,shape=Aq,rate=Bq)
  gauss_m       <- rowMeans(apply(beta_sample,1,function(beta) c(beta-muq,vech(tcrossprod(beta)-sigq-tcrossprod(muq)))))
  invgam_m      <- rowMeans(sapply(sigma2_sample,function(sigma2) c(-log(sigma2)+log(Bq)-digamma(Aq),-1/sigma2+Aq/Bq)))
  gauss_grad    <- 0
  invgam_grad   <- 0
  lb            <- 0
  for (j in 1:S) {
    beta        <- beta_sample[j,]
    sigma2      <- sigma2_sample[j]
    logh_v      <- logh(n,p,y,X,beta,mu0,sig0,sigma2,A,B)
    logq_v      <- logq(p,beta,muq,invsigq,sigma2,Aq,Bq)
    f_v         <- logh_v-logq_v
    lb          <- lb+f_v
    gauss_grad  <- gauss_grad+(c(beta-muq,vech(tcrossprod(beta)-sigq-tcrossprod(muq)))-gauss_m)*f_v
    invgam_grad <- invgam_grad+(c(-log(sigma2)+log(Bq)-digamma(Aq),-1/sigma2+Aq/Bq)-invgam_m)*f_v
  }
  gauss_grad    <- gauss_grad/(S-1)
  invgam_grad   <- invgam_grad/(S-1)
  gauss_lambda  <- (1-rho)*gauss_lambda+rho*c(gauss_invfisherinformation%*%gauss_grad)
  invgam_lambda <- (1-rho)*invgam_lambda+rho*c(solve(invgam_fisherinformation,invgam_grad))
  lambda1       <- gauss_lambda[1:p]
  lambda2       <- gauss_lambda[(p+1):length(gauss_lambda)]
  sigq          <- -0.5*solve(matrix(c(crossprod(D_plus,lambda2)),p,p,byrow=FALSE))
  muq           <- c(sigq%*%lambda1)
  if (any(eigen(sigq)$values<0)) {
    gauss_lambda <- gauss_lambda_old
    sigq         <- sigq_old
    muq          <- muq_old
  } else {
    invsigq      <- solve(sigq)
  }
  Aq            <- invgam_lambda[1]
  Bq            <- invgam_lambda[2]
  if (Aq<0) Aq  <- Aq_old
  if (Bq<0) Bq  <- Bq_old
  invgam_lambda <- c(Aq,Bq)
  lb <- lb/S
  lbc <- c(lbc,lb)
  muq_storage[,i] <- muq
  sigma2_storage[,i] <- invgam_lambda
  cat('lb: ',lb,'\n')
}


# plot
par(mfrow=c(4,2))
plot(lbc[100:length(lbc)],type='l',xlab='# of iterations',ylab='LB',main='Convergence of the Lower Bound')
plot(muq_storage[1,],type='l',xlab='# of iterations',ylab=expression(mu[q1]),main=expression(paste('Convergence of ',mu[q1])))
plot(muq_storage[2,],type='l',xlab='# of iterations',ylab=expression(mu[q2]),main=expression(paste('Convergence of ',mu[q2])))
plot(muq_storage[3,],type='l',xlab='# of iterations',ylab=expression(mu[q3]),main=expression(paste('Convergence of ',mu[q3])))
plot(muq_storage[4,],type='l',xlab='# of iterations',ylab=expression(mu[q4]),main=expression(paste('Convergence of ',mu[q4])))
plot(sigma2_storage[1,],type='l',xlab='# of iterations',ylab=expression(A[q]),main=expression(paste('Convergence of ',A[q])))
plot(sigma2_storage[2,],type='l',xlab='# of iterations',ylab=expression(B[q]),main=expression(paste('Convergence of ',B[q])))



