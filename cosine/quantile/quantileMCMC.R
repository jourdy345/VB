require('ghyp')
require('mvtnorm')
fRho <- function(Yt,x,const1,tau,theta,J) {
	exp(const1*(Yt-x)-sum((exp(Yt*(1:J))-exp(x*(1:J)))*theta^2)/(2*tau)-x+Yt)
}
MCMC <- function(y,x,W,prior,J,p,n_sample,burnin,thinin) {
	n <- length(y)
	if (!is.matrix(W)) W <- as.matrix(W)
	A             <- prior$A
	B             <- prior$B
	muBeta        <- prior$muBeta
	SigmaBeta     <- prior$SigmaBeta
	w0            <- prior$w0
	SigmaBeta_inv <- solve(SigmaBeta)
	SigmaBeta_inv_muBeta <- as.vector(SigmaBeta_inv%*%muBeta)
	if (any(x<0 | x>1)) {
		xmax   <- max(x)
		xmin   <- min(x)
		x_adj  <- (x-xmin)/(xmax-xmin)
		varphi <- sqrt(2/(xmax-xmin))*cos(outer(x_adj,pi*(1:J)))
	} else {
		varphi <- sqrt(2)*cos(outer(x,pi*(1:J)))
	}
	const1   <- J * (J + 1) / 4 - w0
	#----initialize----#
	beta   <- rnorm(ncol(W))
	theta  <- rnorm(J)
	gamma  <- 0.01
	tau    <- 1/rexp(1)
	Aq     <- A+0.5*J
	nup    <- (1-2*p)/(p*(1-p))
	taup2  <- 2/(p*(1-p))
	myRgig <- function(chi) ghyp::rgig(1,0.5,chi,0.25*taup2)
	cat("Burning in...\n")
  count <- 0
	for (t in 1:burnin) {
		count  <- count+1
		cat(count, " iterations (Burn-in)\n")
		Estar  <- exp((1:J)*gamma)
		E      <- Estar/tau
		#----update u----#
		u      <- sapply((y-as.vector(W%*%beta)-as.vector(varphi%*%theta))^2/taup2,myRgig)
		#----update beta----#
		sig.t  <- solve(crossprod(W*sqrt(1/u))/taup2+SigmaBeta_inv)
		mu.t   <- sig.t%*%(colSums(W*(y-varphi%*%theta-nup*u)/u)/taup2+SigmaBeta_inv_muBeta)
		beta   <- as.vector(mvtnorm::rmvnorm(1,mean=mu.t,sigma=sig.t))
		#----update tau----#
		tau    <- 1/rgamma(1,Aq,rate=0.5*sum(Estar*theta^2)+B)
		#----update theta----#
		sig.t  <- solve(crossprod(varphi*sqrt(1/u))/taup2+diag(E))
		mu.t   <- sig.t%*%colSums(varphi*((y-as.vector(W%*%beta)-nup*u)/u))/taup2
		theta  <- as.vector(mvtnorm::rmvnorm(1,mean=mu.t,sigma=sig.t))
		gcan   <- rexp(1)
		rho    <- fRho(gcan,gamma,const1,tau,theta,J)
		gamma  <- gamma+(gcan-gamma)*(runif(1)<rho)
	}

	uList     <- list(u)
	thetaList <- list(theta)
	betaList  <- list(beta)
	tauList   <- c(tau)
	gammaList <- c(gamma)

	count <- 1
	for (j in 2:n_sample) {
		count <- count + 1
		cat(count, " samples collected...\n")
		for (t in 1:thinin) {
			#----update u----#
			u      <- sapply((y-as.vector(W%*%beta)-as.vector(varphi%*%theta))^2/taup2,myRgig)
			#----update beta----#
			Estar  <- exp((1:J)*gamma)
			E      <- Estar/tau
			sig.t  <- solve(crossprod(W*sqrt(1/u))/taup2+SigmaBeta_inv)
			mu.t   <- sig.t%*%(colSums(W*(y-varphi%*%theta-nup*u)/u)/taup2+SigmaBeta_inv_muBeta)
			beta   <- as.vector(rmvnorm(1,mean=mu.t,sigma=sig.t))
			#----update theta----#
			sig.t  <- solve(crossprod(varphi*sqrt(1/u))/taup2+diag(E))
			mu.t   <- sig.t%*%colSums(varphi*((y-as.vector(W%*%beta)-nup*u)/u))/taup2
			theta  <- as.vector(rmvnorm(1,mean=mu.t,sigma=sig.t))
			#----update tau----#
			tau    <- 1/rgamma(1,Aq,rate=0.5*sum(Estar*theta^2)+B)
			gcan   <- rexp(1)
			rho    <- fRho(gcan,gamma,const1,tau,theta,J)
			gamma  <- gamma+(gcan-gamma)*(runif(1)<rho)
		}
		uList[[j]] <- u
		thetaList[[j]] <- theta
		betaList[[j]] <- beta
		tauList <- c(tauList, tau)
		gammaList <- c(gammaList, gamma)
	}
	u     <- do.call('cbind', uList)
	theta <- do.call('cbind', thetaList)
	beta  <- do.call('cbind', betaList)
	list(u=u,theta=theta,beta=beta,gamma=gammaList,tau=tauList,varphi=varphi)
}

x         <- runif(1000)
y         <- sin(2*pi*x)-log(x)
A         <- 0.01
B         <- 0.01
muBeta    <- 0
SigmaBeta <- matrix(1)
w0        <- 1
W         <- matrix(1,nr=length(y),nc=1)
priors    <- list(A=A,B=B,muBeta=muBeta,SigmaBeta=SigmaBeta,w0=w0)
res       <- MCMC(y,x,W,priors,30,0.5,1000,5000,30)
# fit       <- res$varphi%*%res$