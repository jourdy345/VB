library(mvtnorm)
library(MCMCpack)
logh <- function(n,alpha,b,y,lam,phi,psi,sig,sig.a2,sig.l2,sig.p2) {
  v <- -n*lam/2-sig/2
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


# initialize hyperparameters
sig.a2 <- 10
sig.l2 <- 10
sig.p2 <- 10


# initialize variational parameters
mu.q            <- rep(0,d)
# sig.q           <- crossprod(matrix(runif(d*d),d,d))
sig.q           <- rwish(204,diag(d))
invsig.q        <- chol2inv(chol(sig.q))
T_mat_transpose <- chol(invsig.q)
T_mat           <- t(T_mat_transpose)
T_mat_prime     <- T_mat
diag(T_mat_prime) <- log(diag(T_mat))
# S <- 200
nIter        <- 10000 # number of iterations
S            <- 50
count        <- 0
lbc          <- c()
# mu.q.storage <- matrix(0,n+3,niter)
mu.q.storage <- list()
check.converge <- c()
lb.avg <- c()
while (TRUE) {
  # mu.q.old <- mu.q
  # T_mat_old <- T_mat
  # T_mat_transpose_old <- T_mat_transpose
  # T_mat_prime_old <- T_mat_prime



  count <- count+1
  rho   <- 1/(5+count)
  eps   <- 1.0e-6
  # if (count==1) s <- rnorm(n+3,0,3)
  # else s     <- rnorm(n+3)
  s.samples     <- matrix(rnorm((n+3)*S),n+3,S)
  g.mu <- 0
  g.T  <- 0
  for (j in 1:S) {
    s     <- s.samples[,j]
    theta <- mu.q+solve(T_mat_transpose,s)
    b     <- theta[1:n]
    alpha <- theta[n+1]
    lam   <- theta[n+2]
    psi   <- theta[n+3]
    phi   <- exp(psi)/(1+exp(psi))
    sig   <- exp(alpha)
    Ts    <- c(T_mat%*%s)
    g.b1  <- -(1-phi*phi)*b[1]+phi*(b[2]-phi*b[1])-sig/2+sig/2*y[1]*y[1]*exp(-lam-sig*b[1])
    g.b   <- phi*(b[3:n]-phi*b[2:(n-1)])-(b[2:(n-1)]-phi*b[1:(n-2)])-sig/2+sig/2*y[2:(n-1)]*y[2:(n-1)]*exp(-lam-sig*b[2:(n-1)])
    g.bn  <- -(b[n]-phi*b[n-1])-sig/2+sig/2*y[n]*y[n]*exp(-lam-sig*b[n])
    g.a   <- 0.5*sum(y*y*b*exp(alpha-lam-sig*b))-sig/2*sum(b)-alpha/sig.a2
    g.l   <- -n/2+0.5*sum(y*y*exp(-lam-sig*b))-lam/sig.l2
    g.p   <- exp(psi)/((1+exp(psi))*(1+exp(psi)))*(phi*b[1]*b[1]-phi/(1-phi*phi)+sum((b[2:n]-phi*b[1:(n-1)])*b[1:(n-1)]))-psi/sig.p2
    g.mu  <- g.mu+c(g.b1,g.b,g.bn,g.a,g.l,g.p)+Ts
  }
  g.mu <- g.mu/S

  # update mu
  mu.q  <- mu.q+rho*g.mu

  # update T_mat
  g.T   <- -tcrossprod(solve(T_mat_transpose,s),solve(T_mat,g.mu))
  diag(g.T) <- diag(g.T)*diag(T_mat)
  T_mat_prime <- T_mat_prime+rho*g.T
  T_mat <- T_mat_prime
  diag(T_mat) <- exp(diag(T_mat_prime))
  T_mat <- T_mat*lower.tri(T_mat)
  T_mat_transpose <- t(T_mat)

  ################# Comment the below portion to revert ########################
  # sig.q <- tcrossprod(T_mat)
  # if (any(is.infinite(sig.q)) || any(eigen(sig.q)$values < 0)) {
  #   T_mat <- T_mat_old
  #   T_mat_transpose <- T_mat_transpose_old
  #   T_mat_prime <- T_mat_prime_old
  #   mu.q <- mu.q.old
  # }







  lb <- logh(n,alpha,b,y,lam,phi,psi,sig.a2,sig.l2,sig.p2)+(n+3)/2*log(2*pi)-determinant(T_mat)$modulus[1]+0.5*sum(s*s)
  lbc <- c(lbc,lb)
  mu.q.storage[[count]] <- mu.q
  cat('count: ',count,', lb: ',lb,'\n')
  if (count %% 10 == 0) {
    lb_bar <- mean(lbc[(count-10):count])
    lb.avg <- c(lb.avg,lb_bar)
    lb.max <- max(lb.avg)
    if (lb_bar < lb.max) {
      check.converge <- c(check.converge,'below')
    }
    if ((all(check.converge == 'below')) && (length(check.converge)>6)) {
      break
    } else {
      check.converge <- c()
    }
  }
}



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
