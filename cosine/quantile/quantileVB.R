#--------------------------------------------------------------------------------------#
#--Generalized Inverse Gaussian distribution-------------------------------------------#
#--Usage: dgig(x, lambda = 1, chi = 1, psi = 1, logvalue=FALSE)                        #
#-------- pgig(q, lambda = 1, chi = 1, psi = 1, ...)                                   #
#-------- rgig(n, lambda = 1, chi = 1, psi = 1)                                        #
#-------- Egig(lambda, chi, psi, func = c('x', 'logx', '1/x', 'var'), check.pars=TRUE) #
#--Note:  Egig returns the expected value of each functional form ---------------------#
#---------aq2 = chi, bq2 = psi --------------------------------------------------------#
#--------------------------------------------------------------------------------------#
require('ghyp')


#---Expectation of folded-normal---#
efnorm <- function(mu,sigma,sigma2) {
  sigma*sqrt(2/pi)*exp(-mu^2/(2*sigma2))+mu*(1-2*pnorm(-mu/sigma))
}

#-------MGF of folded-normal-------#
mfnorm <- function(mu,sigma,sigma2,t) {
  exp(sigma2*j^2/2+mu*t)*(1-pnorm(-mu/sigma-sigma*t)) + exp(sigma2*j^2/2-mu*t)*(1-pnorm(mu/sigma-sigma*t))
}

#-----------Lower bound------------#
lowerBound <- function(y,W,varphi,quant,J,A,B,Aq,Bq,aq2,bq2,muBeta,SigmaBeta_inv,muBetaq,SigmaBetaq,muThetaq,SigmaThetaq) {
  n <- length(y)
  p <- dim(W)[2]
  taup2 <- 2/(quant*(1-quant))
  nup   <- (1-2*quant)/(quant*(1-quant))
  t1    <- SigmaBeta_inv %*% SigmaBetaq


  res <- -n/2*log(2*pi*taup2) - 0.5*sum(1/taup2*Egig(rep(0.5,n),aq2,rep(bq2,n),func='1/x')*(y-W%*%muBetaq-varphi%*%muThetaq)^2+diag(W%*%tcrossprod(SigmaBetaq,W))+diag(varphi%*%tcrossprod(SigmaThetaq,varphi)))
  -Egig(rep(0.5,n),aq2,rep(bq2,n),func='x') - p/2*log(2*pi)+0.5*determinant(t1)$modulus[1]-0.5*sum((muBetaq-muBeta)*(SigmaBeta_inv%*%(muBetaq-muBeta)))
  -J/2*(log(2*pi)+log(Bq)-psi(Aq))+0.25*J*(J+1)*efnorm(mug,sigmag,sigmag2)+Aq/Bq*sum(mfnorm(mug,sigmag,sigmag2,1:J)*muThetaq^2)
  +log(w0/2)-w0*efnorm(mug,sigmag,sigmag2)+A*log(B)-lgamma(A)-(A+1)*(log(Bq)-digamma(Aq))-B*Aq/Bq
  +0.5*log(sqrt(aq2/bq2))+log(2*besselK(sqrt(aq2*bq2),0.5, expon.scaled=FALSE))+0.5*(aq2*Egig(rep(0.5,n),aq2,rep(bq2,n),func='1/x')+bq2*Egig(rep(0.5,n),aq2,rep(bq2,n),func='x'))
  +(p+J+1)/2*(1+log(2*pi))+0.5*determinant(SigmaThetaq)$modulus[1]+0.5*log(sigmag2)+Aq+log(Bq)+lgamma(Aq)-(1-Aq)*digamma(Aq)

  res
}

