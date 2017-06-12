library(mvtnorm)
library(MCMCpack)
# y <- read.table('~/Desktop/Github/VB/SVI/Code/Gaussian\ VA/VA_svm/exchange.txt')[[1]]

logh <- function(n,alpha,b,y,lam,phi,psi,sig,sig.a2,sig.l2,sig.p2) {
  v <- -n*lam/2-sig/2*sum(b)-0.5*sum(y*y*exp(-lam-sig*b))-
       0.5*sum((b[2:n]-phi*b[1:(n-1)])^2)+0.5*log(1-phi^2)-
       0.5*(1-phi^2)*b[1]^2-alpha^2/(2*sig.a2)-lam^2/(2*sig.l2)-psi^2/(2*sig.p2)
  v     
}

n <- 60
d <- n+3
phi_true <- 0.3
lam_true <- 1.24
sig_true <- -2.73
b1_true <- rnorm(1,0,sqrt(1/(1-phi_true*phi_true)))
b_true <- rep(0,n)
b_true[1] <- b1_true
y <- rep(0,n)
y[1] <- rnorm(1,0,exp(0.5*(lam_true+sig_true*b_true[1])))
for (t in 2:n) {
  b_true[t] <- rnorm(1,phi_true*b_true[t-1],1)
  y[t] <- rnorm(1,0,exp(0.5*(lam_true+sig_true*b_true[t])))
}

# n <- length(y)
# d <- n+3
# initialize hyperparameters

sig.a2 <- 10
sig.l2 <- 10
sig.p2 <- 10


# initialize variational parameters
mu.q     <- rep(0,d)
sig.q    <- diag(d)
invsig.q <- chol2inv(chol(sig.q))
L        <- t(chol(invsig.q))
L_t      <- t(L)
L_p      <- L

count <- 0
lbc   <- c()
# Egmu  <- rep(0,d)
# EgT   <- matrix(0,d,d)
# Edmu  <- rep(0,d)
# dmu   <- rep(0,d)
# EdT   <- matrix(0,d,d)
# dT    <- matrix(0,d,d)

b     <- rep(0,n)
alpha <- 0.5
lam   <- 0.5
psi   <- 0.5
phi   <- exp(psi)/(1+exp(psi))
sig   <- exp(alpha)
g.mu   <- rep(0,d)
g.T    <- matrix(0,d,d)

Egmu  <- rep(0,d)
EgT   <- matrix(0,d,d)
Edmu  <- rep(0,d)
EdT   <- matrix(0,d,d)
dmu   <- rep(0,d)
dT    <- matrix(0,d,d)

nIter <- 90000
mu.q.storage <- list()
check.converge <- c()
epsilon <- 1.0e-06
cnt     <- 0

LBavg   <- 0
LBavg.c <- c()
maxLB  <- -Inf
F <- 2500
for (i in 1:nIter) {
  rho <- 1/(5+i)
  s <- rnorm(d)
  Ls <- solve(L_t,s)
  theta <- as.vector(Ls+mu.q)
  b <- theta[1:n]
  alpha <- theta[n+1]
  lam <- theta[n+2]
  psi <- theta[n+3]
  phi <- 1/(1+exp(-psi))
  sig <- exp(alpha)
  g.b1 <- -(1-phi^2)*b[1]+phi*(b[2]-phi*b[1])-sig/2+
          sig/2*y[1]^2*exp(-lam-sig*b[1])
  g.bt <- phi*(b[3:n]-phi*b[2:(n-1)])-(b[2:(n-1)]-phi*b[1:(n-2)])-sig/2+
          sig/2*y[2:(n-1)]^2*exp(-lam-sig*b[2:(n-1)])
  g.bn <- -(b[n]-phi*b[n-1])-sig/2+sig/2*y[n]^2*exp(-lam-sig*b[n])
  g.a  <- 0.5*sum(y*y*b*exp(alpha-lam-sig*b))-
          sig/2*sum(b)-alpha/sig.a2
  g.l  <- -n/2+0.5*sum(y*y*exp(-lam-sig*b))-lam/sig.l2
  g.p  <- (phi*b[1]^2-phi/(1-phi^2)+sum((b[2:n]-phi*b[1:(n-1)])*b[1:(n-1)]))*exp(psi)/(exp(psi)+1)^2-psi/sig.p2
  g.mu <- c(g.b1,g.bt,g.bn,g.a,g.l,g.p)+Ls
  Egmu <- rho*Egmu+(1-rho)*g.mu^2
  dmu  <- sqrt(Edmu+epsilon)/sqrt(Egmu+epsilon)*g.mu
  Edmu <- rho*Edmu+(1-rho)*dmu^2
  mu.q <- mu.q+dmu
  g.T  <- -tcrossprod(Ls,solve(L,g.mu))
  diag(g.T) <- diag(g.T)*diag(L)
  EgT       <- rho*EgT+(1-rho)*g.T^2
  dT        <- sqrt(EdT+epsilon)/sqrt(EgT+epsilon)*g.T
  EdT       <- rho*EdT+(1-rho)*dT^2
  L_p       <- L_p+dT
  L         <- L_p
  L         <- L*lower.tri(L,diag=TRUE)
  diag(L)   <- exp(diag(L_p))


  lb   <- logh(n,alpha,b,y,lam,phi,psi,sig,sig.a2,sig.l2,sig.p2)+d/2*log(2*pi)-determinant(L)$modulus[1]+0.5*sum(s*s)
  LBavg <- lb/F
  if (is.nan(lb)) break
  mu.q.storage[[i]] <- mu.q
  lbc  <- c(lbc,lb)
  if (i %% F == 0) {
    LBavg.c <- c(LBavg.c,LBavg)
    if (LBavg < maxLB) {
      cnt <- cnt+1
    } else if (LBavg > maxLB) {
      cnt <- 0
    }
    maxLB <- max(maxLB,LBavg)
    cat('LBavg: ',LBavg,', maxLB: ',maxLB,', cnt: ',cnt,'\n')
    LBavg   <- 0
  }
  if (cnt > 4) break
}

mu.q.est <- do.call('cbind',mu.q.storage)
















# lb.avg <- c()
# while (TRUE && any(!is.nan(L))) {
#   count  <- count+1
#   rho    <- 1/(5+count)

#   L_old      <- L
#   s      <- rnorm(d)
#   theta  <- c(L%*%s)+mu.q
#   b      <- theta[1:n]
#   alpha  <- theta[n+1]
#   lam    <- theta[n+2]
#   psi    <- theta[n+3]
#   phi    <- exp(psi)/(1+exp(psi))
#   sig    <- exp(alpha)
#   s1     <- exp(alpha)
#   s2     <- b[2:n]-phi*b[1:(n-1)]
#   s3     <- exp(y-lam-s1*b)
#   g.mu[1:n] <- 0.5*s1*(s3-1)
#   g.mu[2:n] <- g.mu[2:n]-s2
#   g.mu[1:(n-1)] <- g.mu[1:(n-1)]+phi*s2
#   g.mu[1]       <- b[1]*(1-phi^2)
#   g.mu[n+1]     <- -0.5*s1*sum(b)+0.5*s1*sum(b*s3)-alpha/sig.a2
#   g.mu[n+2]     <- -0.5*n+0.5*sum(s3)-lam/sig.l2
#   g.mu[n+3]     <- (sum(s2*b[1:(n-1)])-phi/(1-phi^2)+phi*b[1]^2)*exp(psi)/(exp(psi)+1)^2-psi/sig.p2

#   # g.b1   <- -(1-phi*phi)*b[1]+phi*(b[2]-phi*b[1])-sig/2+sig/2*y[1]*y[1]*exp(-lam-sig*b[1])
#   # g.b    <- phi*(b[3:n]-phi*b[2:(n-1)])-(b[2:(n-1)]-phi*b[1:(n-2)])-sig/2+sig/2*y[2:(n-1)]*y[2:(n-1)]*exp(-lam-sig*b[2:(n-1)])
#   # g.bn   <- -(b[n]-phi*b[n-1])-sig/2+sig/2*y[n]*y[n]*exp(-lam-sig*b[n])
#   # g.a    <- 0.5*sum(y*y*b*exp(alpha-lam-sig*b))-sig/2*sum(b)-alpha/sig.a2
#   # g.l    <- -n/2+0.5*sum(y*y*exp(-lam-sig*b))-lam/sig.l2
#   # g.p    <- exp(psi)/((1+exp(psi))*(1+exp(psi)))*(phi*b[1]*b[1]-phi/(1-phi*phi)+sum((b[2:n]-phi*b[1:(n-1)])*b[1:(n-1)]))-psi/sig.p2
#   # g.mu   <- c(g.b1,g.b,g.bn,g.a,g.l,g.p)

#   mu.q   <- mu.q+rho*g.mu
#   L      <- (L+rho*(tcrossprod(g.mu,s)+diag(1/diag(L))))*lower.tri(L,diag=TRUE)
#   logh_v <- logh(n,alpha,b,y,lam,phi,psi,sig.a2,sig.l2,sig.p2)
#   lb     <- logh_v+determinant(L)$modulus[1]+d/2*log(2*pi)+0.5*sum(s*s)
#   lbc    <- c(lbc,lb)
#   cat('count: ',count,', lb: ',lb,', logh: ',logh_v,'\n')
# }


  # Egmu  <- rho*Egmu+(1-rho)*g.mu^2
  # dmu   <- sqrt((Edmu+eps))/sqrt((Egmu+eps))*g.mu
  # Edmu  <- rho*Edmu+(1-rho)*dmu^2
  # mu.q  <- mu.q+dmu

  # update mu
  # mu.q  <- mu.q+rho*g.mu

  # update T_mat
  # g.T         <- -tcrossprod(solveTs,solve(T_mat,g.mu))
  # diag(g.T)   <- diag(g.T)*diag(T_mat)
  # EgT         <- rho*EgT+(1-rho)*g.T^2
  # dT          <- sqrt((EdT+eps))/sqrt((EgT+eps))*g.T
  # EdT         <- rho*EdT+(1-rho)*dT^2
  # T_mat_prime <- T_mat_prime+dT
  # T_mat       <- T_mat_prime
  # diag(T_mat) <- exp(diag(T_mat_prime))
  
  # T_mat_prime <- T_mat_prime+rho*g.T
  # T_mat <- T_mat_prime
  # diag(T_mat) <- exp(diag(T_mat_prime))
  # # T_mat <- T_mat*lower.tri(T_mat)
  # T_mat_transpose <- t(T_mat)


  ################# Comment the below portion to revert ########################
  # sig.q <- tcrossprod(T_mat)
  # if (any(is.infinite(sig.q)) || any(eigen(sig.q)$values < 0)) {
  #   T_mat <- T_mat_old
  #   T_mat_transpose <- T_mat_transpose_old
  #   T_mat_prime <- T_mat_prime_old
  #   mu.q <- mu.q.old
  # }






  # lb  <- logh(n,alpha,b,y,lam,phi,psi,sig.a2,sig.l2,sig.p2)+(n+3)/2*log(2*pi)-determinant(T_mat)$modulus[1]+0.5*sum(s*s)
  # lb  <- logh_v+(n+3)/2*log(2*pi)-determinant(T_mat)$modulus[1]+0.5*sum(s*s)
  # lbc <- c(lbc,lb)
  # mu.q.storage[[count]] <- mu.q
  # if (count %% 10 == 0) {
  #   lb_bar <- mean(lbc[(count-10):count])
  #   lb.avg <- c(lb.avg,lb_bar)
  #   lb.max <- max(lb.avg)
  #   if (lb_bar < lb.max) {
  #     check.converge <- c(check.converge,'below')
  #   }
  #   if ((all(check.converge == 'below')) && (length(check.converge)>6)) {
  #     break
  #   } else {
  #     check.converge <- c()
  #   }
  # }
# }



# par(mfrow=c(3,4))
# plot(lbc[100:length(lbc)],type='l',xlab='',ylab='LB',main='Convergence of lower bound')
# for (k in 1:7) {
#   plot(mu.q.storage[k,],type='l',xlab='',ylab=substitute(paste('convergence path of ',b[index]),list(index=k)),main=substitute(paste('Convergence of ',b[index]),list(index=k)))
# }

# # var_true <- 1/(1-phi_true*phi_true)
# prec <- tcrossprod(T_mat)
# cov_est <- diag(chol2inv(chol(prec)))

# for (k in 2:5) {
#   curve(dnorm(x,mu.q[k],sqrt(cov_est[k])),col='black',lty=2,lwd=2)
#   curve(dnorm(x,phi_true*b_true[k-1],1),col='red',lwd=2,add=TRUE)
# }
