library(mvtnorm)
library(FNN)            # Get k nearest neighbours
# library(lhs)          # Latin hypercube sampling
# library(mlegp)        # fit guassian process using max likelihood
# library(calibrate)    # To label points on scatterplot
# neigh <- attr(knn(train, test, labels, k = 2, algorithm="cover_tree"),"nn.index") 
# this returns the indices of the k nearest neighbours


#-----------#
# Functions #
#-----------#

# function to compute t(x)%*%A^(-1)%*%x (x is a vector)
quadinv <- function(x,A) sum(x*solve(A,x))

# function to compute trace of A
tr <- function(A) sum(diag(as.matrix(A)))

# function to compute tr(A^(-1)%*%B)
# tri <- function(A,B) sum(diag(as.matrix(solve(A,B)))) 

# rescale variables to have range [-1,1], not necessarily with mean zero though often close
# scle <- function(x) { -1+2*(x-min(x))/(max(x)-min(x)) } 


# NMSE: normalized mean square error
# MNLP: Mean negative log probability
# yp: observations for testing
# prmean: predictive mean
# prvar: predictive variance
NMSE <- function(ypr,prmean,y) {mean((ypr-prmean)^2)/mean((ypr-mean(y))^2)}
MNLP <- function(ypr,prmean,prvar) {0.5*mean( (ypr-prmean)^2/prvar+log(prvar)+log(2*pi) )}


#----------------------------------#
# function to construct T from X,S #
#----------------------------------#

Tfunc <- function(X,S) {
  n <- dim(X)[1]
  m <- dim(S)[1]
  d <- dim(S)[2]
  T <- matrix(0,(n*m),d)
  for (i in 1:n) {
    T[((i-1)*m+1):(i*m),] <- S*matrix(rep(X[i,],m),nrow=m,ncol=d,byrow=TRUE)
  }
  list(T=T)
}


#---------------------#
# Mean of Z and Z^T Z #
#---------------------#

MZ <- function(n,m,T,mulq,siglq) {
  mz <- matrix(0,n,(2*m))
  t1 <- T%*%mulq
  t2 <- exp(-0.5*rowSums((T%*%siglq)*T))
  mz[,1:m] <- matrix(cos(t1)*t2,n,m,byrow=TRUE)
  mz[,(m+1):(2*m)] <- matrix(sin(t1)*t2,n,m,byrow=TRUE) 
  return(mz)
}

MZTZ <- function(n,m,Tm,Tp,mulq,siglq) {
  mztz <- matrix(0,(2*m),(2*m))
  T3m <- matrix(Tm%*%mulq,n,m^2,byrow=TRUE)
  T3p <- matrix(Tp%*%mulq,n,m^2,byrow=TRUE)
  T4m <- matrix(exp(-0.5*rowSums((Tm%*%siglq)*Tm)),n,m^2,byrow=TRUE)
  T4p <- matrix(exp(-0.5*rowSums((Tp%*%siglq)*Tp)),n,m^2,byrow=TRUE)
  mztz[1:m,1:m] <- 0.5*matrix(colSums(T4m*cos(T3m)+T4p*cos(T3p)),m,m,byrow=TRUE)
  mztz[1:m,(m+1):(2*m)] <- 0.5*matrix(colSums(-T4m*sin(T3m)+T4p*sin(T3p)),m,m,byrow=TRUE)
  mztz[(m+1):(2*m),1:m] <- t(mztz[1:m,(m+1):(2*m)])
  mztz[(m+1):(2*m),(m+1):(2*m)] <- 0.5*matrix(colSums(T4m*cos(T3m)-T4p*cos(T3p)),m,m,byrow=TRUE)
  return(mztz)
}


#----------------------------------#
# Function to compute log H(p,q,r) #
#----------------------------------#

h <- function(x,p,q,r) { p*log(x)-q*x^2-log(r+x^(-2)) }
h1 <- function(x,p,q,r) { p/x-2*q*x+2/(r*x^3+x) } #first derivative of h
h2 <- function(x,p,q,r) { -p/x^2-2*q-2*(3*r*x^2+1)/(r*x^3+x)^2 } #second derivative of h
hmaxpt <- function(p,q,r) { sqrt((p*r-2*q+sqrt((p*r-2*q)^2+8*q*r*(p+2)))/(4*q*r)) }

logH <- function(p,q,r){
  mu0 <- hmaxpt(p,q,r)
  sig0 <- (-h2(mu0,p,q,r))^(-0.5)
  hmu0 <- h(mu0,p,q,r)
  sig02 <- sig0*sqrt(2)
  lowerlimit <- (-mu0)/sig02
  integrand <- function(u) {exp(h(mu0+u*sig02,p,q,r)-hmu0)}
  b <- 1
  epsilon <- 1
  while (epsilon > 1.0e-5) {
    b <- 2*b
    if (-b > lowerlimit) {
      epsilon <- max(integrand(b),integrand(-b)) 
    } else {
      epsilon <- integrand(b)
    }
  }

  if (-b > lowerlimit) {
    I0 <- integrate(integrand, lower=-b, upper=b)$value
  } else {
    I0 <- integrate(integrand, lower=lowerlimit, upper=b)$value
  }
  hmu0+log(sig02)+log(I0)
}


#--------------------------------------------------------------------#
# Function to compute expectation of sigma^2,gamma^2,sigma and gamma #
#--------------------------------------------------------------------#

Es2 <- function(C,m,As){ exp(logH(2*m-4,C,As^2)-logH(2*m-2,C,As^2))}
Es <- function(C,m,As){ exp(logH(2*m-3,C,As^2)-logH(2*m-2,C,As^2))}




#-------------#
# Lower bound #
#-------------#

LBC <- function(y,X,A,T,Tm,Tp,As,mul0,sigl0,Csq,muaq,sigaq,mulq,siglq,sigb0,mubq,sigbq) {
  n <- length(y)
  d <- dim(X)[2]
  m <- dim(T)[1]/n
  s <- dim(A)[2]
  t1 <- as.matrix(solve(sigl0,siglq))
  t2 <- as.matrix(solve(sigb0,sigbq))
  mz <- MZ(n,m,T,mulq,siglq) 
  mztz <- MZTZ(n, m, Tm, Tp, mulq, siglq)
  lb = -0.5 * (sum(mztz * sigaq) + sum((mztz - crossprod(mz)) * tcrossprod(muaq)) + sum(crossprod(A) * sigbq)) + sum(log((pnorm(mz %*% muaq + A %*% mubq))^y * (1 - pnorm((mz %*% muaq + A %*% mubq)^(1-y))))) + m * log(m) + 0.5 * (determinant(sigaq)$modulus[1] + determinant(t1)$modulus[1] + determinant(t2)$modulus[1]) + m - 0.5 * (tr(t2) + quadinv(mubq, sigb0) + tr(t1) + quadinv((mulq - mul0), sigl0)) + 0.5 * (s+d)  + log(2 * As) - log(pi) + logH(2*m - 2, Csq, As^2)

  # lb <- ( -n/2*log(2*pi)+d/2+m+m*log(m)+0.5*determinant(as.matrix(sigaq))$modulus[1]
  #         +0.5*determinant(t1)$modulus[1]-0.5*quadinv(mulq-mul0,sigl0)-0.5*tr(t1)
  #         +log(2*As)+log(2*Ag)-2*log(pi)+logH(2*m-2,Csq,As^2)+logH(n-2,Cgq,Ag^2)
  #         -0.5*m*exp(logH(2*m,Csq,As^2)-logH(2*m-2,Csq,As^2))*(crossprod(muaq)+tr(sigaq))
  #         +Csq*exp(logH(2*m,Csq,As^2)-logH(2*m-2,Csq,As^2))
  #         -0.5*exp(logH(n,Cgq,Ag^2)-logH(n-2,Cgq,Ag^2))*(crossprod(y-2*AA)
  #         +sum(MZTZ(n,m,Tm,Tp,mulq,siglq)*(tcrossprod(muaq)+sigaq)))
  #         +Cgq*exp(logH(n,Cgq,Ag^2)-logH(n-2,Cgq,Ag^2))  )
  list(lb=lb)
}


#----------------------------------------------------------------#
# Variational Algorithm (half-Cauchy prior with overrelaxation)  #  
# only applicable to small data where n*m*m does not exceed 10^7 #
#----------------------------------------------------------------#

VARC <- function(y,X,Amat,T,As,mul0,sigl0,sigb0, tol=1.0e-4,fac=1.5,fit=NULL,iter=500) {

  n <- length(y)
  d <- dim(X)[2]
  m <- dim(T)[1]/n
  s <- dim(Amat)[2]
  Tm <- matrix(0,(n*m*m),d)
  Tp <- matrix(0,(n*m*m),d)
  for (i in 1:n) {
    for (r in 1:m) {
      Tm[((i-1)*(m^2)+(r-1)*m+1):((i-1)*(m^2)+r*m),] <- matrix(rep(T[((i-1)*m+r),],m),m,d,byrow=TRUE)-T[((i-1)*m+1):(i*m),]
      Tp[((i-1)*(m^2)+(r-1)*m+1):((i-1)*(m^2)+r*m),] <- matrix(rep(T[((i-1)*m+r),],m),m,d,byrow=TRUE)+T[((i-1)*m+1):(i*m),]
    }
  }

  if (is.null(fit)) {

    # Initialize Csq,Cgq #
    vares <- var(y)
    Csq <- (m-1)*vares

    # Initialize mulq,siglq,muaq #
    mulq <- rep(0.5,d)
    siglq <- 0.5*diag(d)
    muaq <- rep(0.5,2*m)

    # initialize mubq, sigbq #
    EqZ <- MZ(n,m,T,mulq,siglq)
    EqZTZ <- MZTZ(n,m,Tm,Tp,mulq,siglq)
    Eqinvs2 <- exp(logH(2*m,Csq,As^2)-logH(2*m-2,Csq,As^2))
    sigbq <- solve(crossprod(Amat) + solve(sigb0))
    mubq <- crossprod(sigbq,crossprod(Amat, y-(EqZ%*%muaq)))
    muystar = EqZ%*%muaq + Amat %*% mubq + rnorm(n)
    # initialize muaq,sigaq #
    sigaq <- solve(m*Eqinvs2*diag(2*m) + EqZTZ)
    muaq <- as.vector(crossprod(sigaq, crossprod(EqZ,(muystar - Amat%*%mubq))))
    lbold <- -10e7
    lbrecord <- NULL
    count <- 0
    a <- 1
    apre <- 1
    
  } else {
    Csq <- fit$Csq
    mulq <- fit$mulq
    siglq <- fit$siglq
    muaq <- fit$muaq
    sigaq <- fit$sigaq
    mubq <- fit$mubq
    sigbq <- fit$sigbq
    lbold <- fit$lb
    lbrecord <- fit$lbrecord
    count <- dim(fit$lbrecord)[1]
    mulqpre <- fit$mulq
    siglq <- fit$siglq
    apre <- fit$apre
  }

  AL <- sigaq + tcrossprod(muaq)
  dif <- 10
  DIFF <- 1

  while ( dif>tol & count<iter ) {

    count <- count+1
    if (count>1) {a <- fac*apre}

    # update mulq,siglq #
    A <- as.vector(AL[1:m,1:m])
    B <- as.vector(AL[1:m,(m+1):(2*m)])
    C <- as.vector(AL[(m+1):(2*m),(m+1):(2*m)])

    T1 <- T%*%mulq
    T2 <- rep(-y+Amat%*%mubq,rep(m,n))*exp(-0.5*rowSums((T%*%siglq)*T))
    F1 <- -crossprod(as.vector(T2*(muaq[1:m]*cos(T1)+muaq[(m+1):(2*m)]*sin(T1)))*T,T)
    F3 <- 2*colSums(as.vector(T2*(muaq[(m+1):(2*m)]*cos(T1)-muaq[1:m]*sin(T1)))*T)

    T3m <- Tm%*%mulq
    T3p <- Tp%*%mulq
    T4m <- exp(-0.5*rowSums((Tm%*%siglq)*Tm))
    T4p <- exp(-0.5*rowSums((Tp%*%siglq)*Tp))

    F2 <- -0.25*( crossprod(as.vector(T4m*((A+C)*cos(T3m)+2*B*sin(T3m)))*Tm,Tm) 
                + crossprod(as.vector(T4p*((A-C)*cos(T3p)+2*B*sin(T3p)))*Tp,Tp) )
    F4 <- 0.5*colSums( as.vector(T4m*(-(A+C)*sin(T3m)+2*B*cos(T3m)))*Tm 
                     + as.vector(T4p*((C-A)*sin(T3p)+2*B*cos(T3p)))*Tp  )

    temp1 <- solve(sigl0)+(F1+F2)
    temp2 <- solve(sigl0,mul0-mulq)-0.5*(F3+F4)

    siglqnew <- solve((1-a)*solve(siglq)+a*temp1)
    EV <- eigen(siglqnew,only.values=TRUE)$values

    while (isSymmetric(siglqnew,tol=1.0e-10)==FALSE | any(Re(EV)<0) ) {
      a <- 2/3*a
      siglqnew <- solve((1-a)*solve(siglq)+a*temp1)
      EV <- eigen(siglqnew,only.values=TRUE)$values
    }

    siglq <- siglqnew
    mulq <- mulq+a*crossprod(siglq,temp2)
    apre <- a
    EqZ <- MZ(n,m,T,mulq,siglq)
    EqZTZ <- MZTZ(n,m,Tm,Tp,mulq,siglq)

    # update muystar
    muystar_temp = EqZ %*% muaq + Amat %*% mubq
    muystar = muystar_temp + dnorm(muystar_temp)/(((pnorm(muystar_temp))^y)*((pnorm(muystar_temp)-1)^(1-y)))
    
    # update muaq,sigaq #
    sigaq <- solve(m * Eqinvs2 * diag(2*m) + EqZTZ)
    muaq <- as.vector(crossprod(sigaq, crossprod(EqZ,(muystar - Amat%*%mubq))))

    # update Csq #
    
    Csq <- m/2*as.numeric(crossprod(muaq)+tr(sigaq))
    lb <- LBC(y,X,Amat,T,Tm,Tp,As,mul0,sigl0,Csq,muaq,sigaq,mulq,siglq,sigb0,mubq,sigbq)$lb

    DIFF <- lb-lbold
    print(DIFF)
    if (count==1 & DIFF<0) {
      cat('DIVERGE',"\n")
    break
    }

    if (DIFF>0){
      mulqpre <- mulq
      siglqpre <- siglq
    } else {
      count <- count+1
      a <- 1
      siglqnew <- solve(temp1)
      EV <- eigen(siglqnew,only.values=TRUE)$values

      while (isSymmetric(siglqnew,tol=1.0e-10)==FALSE | any(Re(EV)<0) ) {
        a <- 2/3*a
        siglqnew <- solve((1-a)*solve(siglq)+a*temp1)
        EV <- eigen(siglqnew,only.values=TRUE)$values
      }

      siglq <- siglqnew
      mulq <- mulqpre+a*crossprod(siglq,temp2)
      apre <- a
      EqZ <- MZ(n,m,T,mulq,siglq)
      EqZTZ <- MZTZ(n,m,Tm,Tp,mulq,siglq)

      # update muystar
      muystar_temp = EqZ %*% muaq + Amat %*% mubq
      muystar = muystar_temp + dnorm(muystar_temp)/(((pnorm(muystar_temp))^y)*((pnorm(muystar_temp)-1)^(1-y)))
    
      # update muaq,sigaq #
      sigaq <- solve(m * Eqinvs2 * diag(2*m) + EqZTZ)
      muaq <- as.vector(crossprod(sigaq, crossprod(EqZ,(muystar - Amat%*%mubq))))

      # update mubq,sigbq #
      sigbq <- solve(crossprod(Amat) + solve(sigbq))
      mubq <- as.vector(crossprod(sigbq, crossprod(Amat, muystar-EqZ%*%muaq)))

      # update Cgq,Csq #
    
      Csq <- m/2*as.numeric(crossprod(muaq)+tr(sigaq))

      lb <- LBC(y,X,Amat,T,Tm,Tp,As,mul0,sigl0,Csq,muaq,sigaq,mulq,siglq,sigb0,mubq,sigbq)$lb
    }

    lbrecord <- rbind(lbrecord,c(count,lb))
    dif <- abs((lb-lbold)/lb)
    lbold <- lb
    cat(count,lb,round(mulq,1),Csq,dif,DIFF,apre,"\n")
  }

  list(Csq=Csq,muaq=muaq,sigaq=sigaq, mubq = mubq, sigbq = sigbq, muystar = muystar,
  mulq=mulq,siglq=siglq,lb=lb,lbrecord=lbrecord,apre=apre)
}


library(mixtools)

#DATA
set.seed(123)
n_=100 #number of i
N_=500 # i * j
id=sample(1:n_,N_,replace=TRUE)
#id=rep(1:n,each=5)
#b=rnorm(n,0,sqrt(2)) # normal random effect with sigma=2
#b=c(rnorm(n/4),rnorm(n/2,0,4),rnorm(n/4,1,2)) # DP random effect

p_=c(.2,.5,.3)
m_=c(-2,0,1.5)
s_=c(.5,1,.5)
b_=rnormmix(n_,p_,m_,s_)
xobs=runif(N_)

f1=function(x)2*x-1
f2=function(x)0.5*exp(x*2)-2
f3=function(x) log(2*x+0.5)
f4 = function(x) tanh(4*x-2)

# plot(xobs,f1(xobs))
# plot(xobs,f2(xobs))
# plot(xobs,f3(xobs))
# plot(xobs,f4(xobs))

y1=rbinom(length(xobs),1,pnorm(f1(xobs)+b[id]))
y2=rbinom(length(xobs),1,pnorm(f2(xobs)+b[id]))
y3=rbinom(length(xobs),1,pnorm(f3(xobs)+b[id]))
y4=rbinom(length(xobs),1,pnorm(f4(xobs)+b[id]))


X = cbind(1, xobs)
Zmat = as.matrix(rep(1, 500))
d = dim(X)[2]
n = dim(X)[1]
m = 10
S = rmvnorm(m, mu = rep(0, d), sigma = diag(d))
T = Tfunc(X = X, S = S)$T
s = 1
fit = VARC(y = y1, X = X, Amat = Zmat, T = T, As = 25, mul0 = rep(0, d), sigl0 = 10 * diag(d), sigb0 = 10 * diag(s), fac = 1.5)


res = sapply(fit$muystar, categorize)


### simulation
sim_GPprobit = function(FUN, intercept = TRUE) {
  if (!require('mixtools')) {
    ans = readline(prompt = "We need to install {mixtools}. Permit installation? (y/n)\n")
    if (ans == 'y') {
      install.packages('mixtools')
      library(mixtools)
    } else {
      stop('Cannot proceed without {mixtools}. Aborting...\n')
    }
  }

  set.seed(123)
  n_=100 #number of i
  N_=500 # i * j
  id=sample(1:n_,N_,replace=TRUE)
  p_=c(.2,.5,.3)
  m_=c(-2,0,1.5)
  s_=c(.5,1,.5)
  b_=rnormmix(n_,p_,m_,s_)
  xobs=runif(N_)
  y = rbinom(length(xobs),1,pnorm(FUN(xobs)+b[id]))
  if (intercept == TRUE) {
    X = cbind(1, xobs)
  } else {
    X = as.matrix(xobs)
  }

  categorize = function(x) {
    if (x < 0) {
      x = 0
    } else {
      x = 1
    }
  }

}
