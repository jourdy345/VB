library(mvtnorm)    # for multivariate Gaussian distribution
library(matrixcalc) # for duplication matrix and elimination matrix
library(MCMCpack)   # for vech and xpnd
library(Matrix)     # for block-diagonal operator
logh <- function(n,alpha,b,y,lam,phi,psi,sig.a2,sig.l2,sig.p2) {
  v <- -0.5*n*lam-0.5*exp(alpha)*sum(b)-0.5*sum(y*y*exp(-lam-exp(alpha)*b))-
       0.5*sum((b[2:length(b)]-b[1:(length(b)-1)])*(b[2:length(b)]-b[1:(length(b)-1)]))-0.5*log(1-phi*phi)-
       0.5*(1-phi*phi)*b[1]*b[1]-alpha*alpha/(2*sig.a2)-lam*lam/(2*sig.l2)-psi*psi/(2*sig.p2)
  v     
}

n <- 200
d <- n+3
phi_true <- 0.3
lam_true <- 1.24
sig_true <- -1.73
b1_true <- rnorm(1,0,sqrt(1/(1-phi_true*phi_true)))
b_true <- rep(0,n)
b_true[1] <- b1_true
y <- rep(0,n)
y[1] <- rnorm(1,0,exp(0.5*(lam_true+sig_true*b_true[1])))
for (t in 2:n) {
  b_true[t] <- rnorm(1,phi_true*b_true[t-1],1)
  y[t] <- rnorm(1,0,exp(0.5*(lam_true+sig_true*b_true[t])))
}


# initialize hyperparameters
sig.a2 <- 1
sig.l2 <- 1
sig.p2 <- 1


# initialize variational parameters
mu.q  <- rep(0,n+3)
sig.q <- crossprod(matrix(runif(d*d),d,d))
invsig.q <- chol2inv(chol(sig.q))

D_mat <- duplication.matrix(d)
D_plus <- ginv(D)
count <- 0
lb <- 0
lbc <- c()
S <- 40
nIter <- 4000

gauss.lambda1 <- c(solve(sig.q,mu.q))
guass.lambda2 <- -0.5*crossprod(D,as.vector(invsig.q))
gauss.lambda <- c(gauss.lambda1,gauss.lambda2)


mu.q.storage <- matrix(0,d,nIter)
count <- 0
# start inference
for (i in 1:nIter) {
  count <- count+1
  cat(count,'th iteration\n')
  gauss.lambda.old <- gauss.lambda
  mu.q.old <- mu.q
  sig.q.old <- sig.q
  rho <- 1/(5+i)
  M_mat <- 2*D_plus%*%(mu.q%x%diag(d))
  S_mat <- 2*D_plus%*%tcrossprod((sig.q%x%sig.q),D_plus)
  aux   <- solve(S_mat,M_mat)
  gauss.invfisherinformation <- rbind(cbind(invsig.q+crossprod(M_mat,aux),-t(aux)),cbind(-aux,solve(S_mat)))

  theta.sample <- rmvnorm(S,mean=mu.q,sigma=sig.q)
  gauss_m      <- rowMeans(apply(theta.sample,1,function(theta) c(theta-mu.q,vech(tcrossprod(theta)-sig.q-tcrossprod(mu.q)))))
  gauss.grad   <- 0
  lb           <- 0
  for (j in 1:S) {
    theta  <- theta.sample[j,]
    b      <- theta[1:n]
    alpha  <- theta[n+1]
    lam    <- theta[n+2]
    psi    <- theta[n+3]
    phi    <- exp(psi)/(1+exp(psi))
    sig    <- exp(alpha)
    logh.v <- logh(n,alpha,b,y,lam,phi,psi,sig.a2,sig.l2,sig.p2)
    logq.v <- dmvnorm(theta,mean=mu.q,sigma=sig.q,log=TRUE)
    f.v    <- logh.v-logq.v
    lb     <- lb+f.v
    gauss.grad <- gauss.grad+(c(theta-mu.q,vech(tcrossprod(theta)-sig.q-tcrossprod(mu.q)))-gauss.m)*f.v
  }
  gauss.grad <- gauss.grad/(S-1)
  gauss.lambda <- (1-rho)*gauss.lambda+rho*c(gauss.invfisherinformation%*%gauss.grad)
  lambda1 <- gauss.lambda[1:d]
  lambda2 <- gauss.lambda[(d+1):length(gauss.lambda)]
  sig.q   <- -0.5*solve(matrix(c(crossprod(D_plus,lambda2)),d,d,byrow=FALSE))
  mu.q    <- c(sig.q%*%lambda1)
  if (any(eigen(sig.q)$values<0)) {
    gauss.lambda <- gauss.lambda.old
    sig.q        <- sig.q.old
    mu.q         <- mu.q.old
  } else {
    invsig.q     <- chol2inv(chol(sig.q))
  }
  lb <- lb/S
  lbc <- c(lbc,lb)
  mu.q.storage[,i] <- mu.q
  cat('lb: ',lb,'\n')
}






