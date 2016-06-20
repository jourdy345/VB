CI_noshape <- function(object,vphi,Z,size=500,quantile=0.95)
{
	#object is the output of the uncontrainted VB algorithm
	#size is the number of samples to be generated from the posterior distribution
	CIbeta <- rmvnorm(size,mean=object$mubeta.q,sigma = object$sigbeta.q)
	CItheta <- rmvnorm(size,mean=object$mutheta.q,sigma = object$sigtheta.q)
	Fullfx = matrix(0, nrow(Z),size)
	FullZB  = matrix(0, nrow(Z),size)	
	for (i in 1:size)
	{
		# ZB <- Z%*%CIbeta[i,]
		FullZB[,i] <- Z%*%CIbeta[i,]
		fx <- (vphi[,1:length(CItheta[i,])]%*%CItheta[i,])
  		Fullfx[,i] = fx
	}
	CIfx <- apply(Fullfx, 1,quantile, probs = c(0.025,0.975))
	CIZB <- apply(FullZB, 1,quantile, probs = c(0.025,0.975))
	object <- list(CIfx=CIfx , CIZB=CIZB)
}

CI_shape <- function(object,Z,size=500,quantitle=0.95,delta)
{	
	#object is the output of the shape restricted VB algorithm
	#size is the number of samples to be generated from the posterior distribution
	CIbeta <- rmvnorm(size,mean=object$mubeta.q,sigma = object$sigbeta.q)
	CItheta <- rmvnorm(size,mean=object$mutheta.q,sigma = object$sigtheta.q)
	FullIntervalfx = matrix(0,length(object$x.grid),size)
	FullIntervalZB  = matrix(0,length(object$x.grid),size)
	for (i in 1:size)
	{
		firstterm <- (Z)%*%CIbeta[i,]
		FullIntervalZB[,i] <- firstterm
		secondterm <- apply(object$dmats,1,fitfun,CItheta[i,]) 
  		FullIntervalfx[,i] = delta*(secondterm)
	}
	CIfx <- apply(FullIntervalfx, 1,quantile, probs = c(0.025,0.975))
	CIZB <- apply(FullIntervalZB, 1,quantile, probs = c(0.025,0.975))
	object <- list(CIfx=CIfx , CIZB=CIZB)
}

#Alternative way of producing credible interval
CI_shape2 <- function(object)
{
  matfun6<-function(psi,sigma,mu) {
    return(2*sum(diag(psi%*%sigma%*%psi%*%sigma))+4*sum(mu*(psi%*%sigma%*%psi%*%mu)))
  }
varterm = apply(object$dmats,1,matfun6,test$sigtheta.q,object$mutheta.q)
meanterm = delta_con*(term1 + term2)

size = 500
FullInterval = matrix(0,length(test$x.grid),size)
for (i in 1:length(test$x.grid))
{
FullInterval[i,] = rnorm(size,meanterm[i],sqrt(varterm[i]))
}
CI = apply(FullInterval, 1,quantile, probs = c(0.025,0.975))


}