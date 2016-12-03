library(mvtnorm)    # for multivariate Gaussian distribution
library(matrixcalc) # for duplication matrix and elimination matrix
library(MCMCpack)   # for vech and xpnd
library(Matrix)     # for block-diagonal operator
library(R.matlab)

# q(beta) ~ N(muq,sigq)
# q(Om)   ~ W(k_q,V_q)

# log-joint density (log(likelihood*prior))
logh <- function(T,p,Y,F,Om,beta,mu0,invsig0,k,invV) {
  logdetOm <- determinant(Om)$modulus[1]

  v <- 0.5*T*logdetOm-0.5*sum(Om*tcrossprod(Y-apply(F,3,function(X) crossprod(X,beta))))-0.5*T*p*log(2*pi)+
       0.5*determinant(invsig0)$modulus[1]-0.5*sum((beta-mu0)*(invsig0)%*%(beta-mu0))+
       (k-p-1)/2*logdetOm-0.5*sum(invV*Om)+k/2*determinant(invV)$modulus[1]-k*p/2*log(2)-
       sum(sapply(1:p,function(x) lgamma((k+1-x)/2)))
  v
}

# log-variational posterior
logq <- function(beta,Om,muq,invsigq,k_q,invV_q,p) {
  v <- 0.5*determinant(invsigq)$modulus[1]-0.5*sum((beta-muq)*(invsigq%*%(beta-muq)))+
       (k_q-p-1)/2*determinant(Om)$modulus[1]-0.5*sum(invV_q*Om)+k_q/2*determinant(invV_q)$modulus[1]-k_q*p/2*log(2)-
       sum(sapply(1:p,function(x) lgamma((k_q+1-x)/2)))
  v
}


data = readMat('Simu_p6.mat')
adj_save = data$adj.save
b0 = data$b0
D0 = data$D0
F = data$F
Fnp = data$Fnp
Fpn = data$Fpn
logm1_save = data$logm1.save
logm2_save = data$logm2.save
logmarg_save = data$logmarg.save
logpost_save = data$logpost.save 
m0 = data$m0
nu0 = data$nu0
nu1 = data$nu1
theta_true = data$theta.true
Y = data$Y
z1_save = data$z1.save

n = dim(F)[1] ; p = dim(F)[2] ; T = dim(F)[3]

n1 = 6    # True regression coefficients
beta1 = c(1.3, 0, -0.5, rep(0,n1))   # y_{1} = 1.3x_{1} - 0.5x_{3} + e_{1}
beta2 = c(0.9, -.3,  0.5, rep(0,n1))  # y_{2} = 0.9z_{1} - 0.3x_{2} + 0.5x_{3} + e_{2}
beta3 = c(1, 0.5,  0.7, rep(0,n1))    # y_{3} = x_{1} + 0.5x_{2} + 0.7x_{3} + e_{3}
beta4 = c(0.8, -0.6, rep(0,n1))       # y_{4} = 0.8x_{4} - 0.6x_{5} + e_{4}
beta5 = c(1, 0.7, rep(0,n1))           # y_{5} = x_{4} + 0.7x_{5} + e_{5}
beta6 = c(1.1, 0.6, rep(0,n1))        # y_{6} = 1.1x_{4} - 0.6x_{5} + e_{6}
beta_true = c(beta1, beta2, beta3, beta4,beta5, beta6)
z1_true = which(abs(beta_true) > 0.001)  # z1_true: true subset of variables;

phi = .6   # True error covariance matrix
Cov_true = 1 / (1 - phi^2) * toeplitz(phi^(0:(p-1)))  
Omega_true = solve(Cov_true) ; adj_true = (abs(Omega_true)>.001)*1  # adj_true: true adjacency matrix;
# Omega = precision matrix, V^{-1}

m0 = matrix(0,n,1)
C0 = 100 * diag(n)
b0 = 3 ; D0 = 0.0001 * diag(p) ; nu0 = 0.01 * matrix(1,n,1) ; nu1= 10 * matrix(1,n,1)
# nu0 = tau_0, nu1 = tau_1

Fpn = NULL
Fnp = NULL

for(t in 1:T){
  Fpn = cbind(Fpn, t(F[,,t]))
  Fnp = cbind(Fnp,F[,,t])
}

# initialize hyperparameters
k <- 3
V <- crossprod(matrix(runif(p*p),p,p))
invV <- chol2inv(chol(V))
mu0 <- rep(0,n)
sig0 <- diag(n)
invsig0 <- chol2inv(chol(sig0))

# initialize variational parameters
k_q <- 10
V_q <- V
invV_q <- chol2inv(chol(V_q))
muq <- rep(0,n)
sigq <- diag(n)
invsigq <- chol2inv(chol(sigq))

D_mat      <- duplication.matrix(n)
D_plus     <- ginv(D_mat)
D_mat_wish <- duplication.matrix(p)

count <- 0
niter <- 2000 # number of iterations
S <- 300 # number of Monte-Carlo samples
lbc <- c()

# natural parameters for beta
gauss_lambda1 <- c(solve(sigq,muq))
guass_lambda2 <- -0.5*crossprod(D_mat,as.vector(invsigq))
gauss_lambda  <- c(gauss_lambda1,guass_lambda2)

wish_lambda <- c(vech(V_q),k_q)

for (i in niter) {
  count <- count+1
  gauss_lambda_old <- gauss_lambda
  muq_old           <- muq
  sigq_old          <- sigq
  wish_lambda_old  <- wish_lambda
  rho   <- 1/(5+i)
  M_mat <- 2*D_plus%*%(muq%x%diag(n))
  S_mat <- 2*D_plus%*%tcrossprod((sigq%x%sigq),D_plus)
  aux   <- solve(S_mat,M_mat)
  gauss_invfisherinformation <- rbind(cbind(invsigq+crossprod(M_mat,aux),-t(aux)),cbind(-aux,solve(S_mat)))
  wish_auxmat <- 0.5*crossprod(D_mat_wish,invV_q%x%invV_q)%*%D_mat_wish
  

  beta_samples <- rmvnorm(S,mean=muq,sigma=sigq)
  Om_samples   <- sapply(1:S,function(x) rwish(k_q,V_q),simplify='array')
  gauss_m      <- rowMeans(apply(beta_samples,1,function(beta) c(beta-muq,vech(tcrossprod(beta)-sigq-tcrossprod(muq)))))
  wish_m       <- rowMeans(apply(Om_samples,3,function(Om) wish_auxmat%*%vech(Om)-k_q/2*vech(invV_q)))
  gauss_grad   <- 0
  wish_grad    <- 0
  lb           <- 0
  
  for (j in 1:S) {
    beta       <- beta_samples[j,]
    Om         <- Om_samples[,,j]
    logh_v     <- logh(T,p,Y,F,Om,beta,mu0,invsig0,k,invV)
    logq_v     <- logq(beta,Om,muq,invsigq,k_q,invV_q,p)
    f_v        <- logh_v-logq_v
    lb         <- lb+f_v
    gauss_grad <- gauss_grad+(c(beta-muq,vech(tcrossprod(beta)-sigq-tcrossprod(muq)))-gauss_m)*f_v
    wish_grad  <- wish_grad+(as.vector(wish_auxmat%*%vech(Om))-k_q/2*vech(invV_q)-wish_m)*f_v
  }
  
  gauss_grad   <- gauss_grad/(S-1)
  wish_grad    <- wish_grad/(S-1)
  gauss_lambda <- (1-rho)*gauss_lambda+rho*c(gauss_invfisherinformation%*%gauss_grad)
  wish_lambda  <- (1-rho)*wish_lambda+rho*wish_grad
  lambda1      <- gauss_lambda[1:n]
  lambda2      <- gauss_lambda[(n+1):length(gauss_lambda)]
  sigq         <- -0.5*solve(matrix(c(crossprod(D_plus,lambda2)),n,n,byrow=FALSE))
  muq          <- c(sigq%*%lambda1)
  if (any(eigen(sigq)$values<0)) {
    gauss_lambda <- gauss_lambda_old
    sigq         <- sigq_old
    muq          <- muq_old
  } else {
    invsigq     <- solve(sigq)
  }
  k_q           <- wish_lambda[length(wish_lambda)]
  V_q           <- xpnd(wish_lambda[1:(length(wish_lambda)-1)],p)
  if (any(eigen(V_q)$values<0)) {
    wish_lambda <- wish_lambda_old
    k_q         <- wish_lambda[length(wish_lambda)]
    V_q         <- xpnd(wish_lambda[1:(length(wish_lambda)-1)])
  } else {
    invV_q <- chol2inv(chol(V_q))
  }
  lb <- lb/S
  lbc <- c(lbc,lb)
  cat('lb: ',lb,'\n')
}




