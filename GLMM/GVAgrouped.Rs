library(statmod)
source("CalculateB.Rs")

###############################################################################

tr <- function(mA) {
   return(sum(diag(mA)))
}

###############################################################################

GVA.grouped.f <- function(vtheta,family,vy,mX,lmZ,gh)  
{
    n <- length(vy)
    p <- ncol(mX)
    m <- length(lmZ) 
    K <- ncol(lmZ[[1]])
    
    Dinds <- K*((1:K)-1) + (1:K) 	
    Rinds <- which(upper.tri(diag(0,K),diag=TRUE))	 
    Linds <- which(lower.tri(diag(0,K),diag=TRUE))	  
    q <- length(Rinds)
    
    ###########################################################################
    
    vbeta <- vtheta[1:p]
    vr    <- vtheta[p+(1:q)]
    mR <- matrix(0,K,K)
    mR[Rinds] <- vr
    mR[Dinds] <- exp(mR[Dinds])
    mSigma.inv <- t(mR)%*%mR

    sl <- 0
    
    vmu.X   <- mX%*%vbeta
    vmu.Z   <- c()
    vsigma2 <- c()    
    mM  <- matrix(0,K,K)
    mRi <- matrix(0,K,K)
    for (i in 1:m) {
        vmui <- vtheta[p+q+(i-1)*(K+q)+(1:K)]
        vri  <- vtheta[p+q+(i-1)*(K+q)+K+(1:q)]
        mRi[Rinds] <- vri
        mRi[Dinds] <- exp(mRi[Dinds])
        mM <- mM + vmui%*%t(vmui) + t(mRi)%*%mRi
        sl <- sl + sum(log(mRi[Dinds]))
        vmu.Z   <- c(vmu.Z,lmZ[[i]]%*%vmui)
        vsigma2 <- c(vsigma2,((lmZ[[i]]%*%t(mRi))^2)%*%matrix(1,K,1))
    }
    veta <- vmu.X + vmu.Z
    
    ###########################################################################    
    
    vB0 <- B0.fun(family,veta,vsigma2,gh) 
    sl <- sl + sum(vy*veta - vB0) 
    sl <- sl + m*sum(log(mR[Dinds])) + 0.5*m*K - 0.5*tr(mM%*%mSigma.inv)
    if (family=="POISSON") {
        sl <- sl - sum(lgamma(vy+1))
    }
    
    return(sl)
}

###############################################################################

GVA.grouped.vg <- function(vtheta,family,vy,mX,lmZ,gh) 
{ 
    n <- length(vy)
    p <- ncol(mX)
    m <- length(lmZ)
    K <- ncol(lmZ[[1]])
    
    Dinds <- K*((1:K)-1) + (1:K) 	
    Rinds <- which(upper.tri(diag(0,K),diag=TRUE))	  
    Linds <- which(lower.tri(diag(0,K),diag=TRUE))	  
    q <- length(Rinds)
    
    ###########################################################################
    
    vbeta <- vtheta[1:p]
    vr    <- vtheta[p+(1:q)]
    mR <- matrix(0,K,K)
    mR[Rinds] <- vr
    mR[Dinds] <- exp(mR[Dinds])
    mSigma.inv <- t(mR)%*%mR    
  
    ###########################################################################
    
    voneK <- matrix(1,K,1)
    
    vmu.X   <- mX%*%vbeta
    vmu.Z   <- c()
    vsigma2 <- c()    
    mM  <- matrix(0,K,K)
    mRi <- matrix(0,K,K)
    for (i in 1:m) {
        vmui <- vtheta[p+q+(i-1)*(K+q)+(1:K)]
        vri  <- vtheta[p+q+(i-1)*(K+q)+K+(1:q)]
        mRi[Rinds] <- vri
        mRi[Dinds] <- exp(mRi[Dinds])
        mM <- mM + vmui%*%t(vmui) + t(mRi)%*%mRi
        vmu.Z   <- c(vmu.Z,lmZ[[i]]%*%vmui)
        vsigma2 <- c(vsigma2,((lmZ[[i]]%*%t(mRi))^2)%*%voneK)
    }
    
    veta <- vmu.X + vmu.Z    
    
    ###########################################################################    
    
    lvB  <- B1234.fun(family,veta,vsigma2)
    vB1 <- lvB[[1]]
    vB2 <- lvB[[2]]; dim(vB2) <- NULL
    vB3 <- lvB[[3]]; dim(vB3) <- NULL
    vB4 <- lvB[[4]]; dim(vB4) <- NULL
    veps <- vy - vB1 
    
    ###########################################################################
    
    vg.vbeta <- t(mX)%*%veps
    mA <- -mR%*%mM
    mA[Dinds] <- m + mA[Dinds]*mR[Dinds]
    vg.mR <- mA[Rinds]
    vg.vtheta <- c(vg.vbeta,vg.mR)
      
    ###########################################################################
    
    vg <- vg.vtheta
    
    count <- 0
    linds <- list()
    lvg.vxi <- list()
    for (i in 1:m) {
        mZi  <- lmZ[[i]]    
        inds <- count + (1:nrow(mZi))
        count <- count + nrow(mZi)
        linds[[i]] <- inds
        vmui <- vtheta[p+q+(i-1)*(K+q)+(1:K)]
        vri  <- vtheta[p+q+(i-1)*(K+q)+K+(1:q)]
        mRi[Rinds] <- vri
        mRi[Dinds] <- exp(mRi[Dinds])
        vB2i <- vB2[inds]
        vepsi <- veps[inds] 
        vg.vmui <- t(mZi)%*%vepsi - mSigma.inv%*%vmui  
        mA <-   -mRi%*%(t(mZi*vB2i)%*%mZi + mSigma.inv)
        mA[Dinds] <- 1 + mA[Dinds]*mRi[Dinds]
        vg.mRi <- mA[Rinds]
        lvg.vxi[[i]] <- matrix(c(vg.vmui,vg.mRi))
        
        vg <- c(vg,lvg.vxi[[i]])
    }
    
    res.vg <- list(vg.0=vg.vtheta,lvg.i=lvg.vxi,vg=vg)
    
    return(res.vg)
}

###############################################################################

GVA.grouped.vgANDmH <- function(vtheta,family,vy,mX,lmZ,gh) 
{ 
    n <- length(vy)
    p <- ncol(mX)
    m <- length(lmZ)
    K <- ncol(lmZ[[1]])
    
    Dinds <- K*((1:K)-1) + (1:K) 	
    Rinds <- which(upper.tri(diag(0,K),diag=TRUE))	  
    Linds <- which(lower.tri(diag(0,K),diag=TRUE))	  
    q <- length(Rinds)
    RDinds <- rep(FALSE,length(Rinds))
    for (i in 1:length(Dinds)) {
        RDinds[i] <- which(Rinds==Dinds[i])
    }
    Rinds2D <- which(upper.tri(diag(0,K),diag=TRUE),arr.ind=T)
    Linds2D <- which(lower.tri(diag(0,K),diag=TRUE),arr.ind=T)
        
    ###########################################################################
    
    vbeta <- vtheta[1:p]
    
    vr    <- vtheta[p+(1:q)]
    mR <- matrix(0,K,K)
    mR[Rinds] <- vr
    mR[Dinds] <- exp(mR[Dinds])
    mSigma.inv <- t(mR)%*%mR    
  
    ###########################################################################
    
    vmu.X   <- mX%*%vbeta
    vmu.Z   <- c()
    vsigma2 <- c()    
    mM  <- matrix(0,K,K)
    mRi <- matrix(0,K,K)
    for (i in 1:m) {
        vmui <- vtheta[p+q+(i-1)*(K+q)+(1:K)]
        vri  <- vtheta[p+q+(i-1)*(K+q)+K+(1:q)]
        mRi[Rinds] <- vri
        mRi[Dinds] <- exp(mRi[Dinds])
        mM <- mM + vmui%*%t(vmui) + t(mRi)%*%mRi
        mZi <- lmZ[[i]]
        vmu.Z   <- c(vmu.Z,mZi%*%vmui)
        vsigma2 <- c(vsigma2,((mZi%*%t(mRi))^2)%*%matrix(1,K,1))
    }
    
    veta <- vmu.X + vmu.Z    
    
    ###########################################################################    
    
    lvB  <- B1234.fun(family,veta,vsigma2)
    vB1 <- lvB[[1]]
    vB2 <- lvB[[2]]; dim(vB2) <- NULL
    vB3 <- lvB[[3]]; dim(vB3) <- NULL
    vB4 <- lvB[[4]]; dim(vB4) <- NULL
    veps <- vy - vB1 
    
    ###########################################################################
        
    res.vg <- GVA.grouped.vg(vtheta,family,vy,mX,lmZ,gh) 
    vg.vtheta <- res.vg$vg.vtheta
    lvg.vxi   <- res.vg$lvg.vxi 
        
    ###########################################################################
    
    mA <- mR%*%mM
    mH.mR.mR <-  matrix(0,q,q)
    for (u in 1:q) {
        i <- Rinds2D[u,1]
        j <- Rinds2D[u,2]
        for (v in 1:q) {
		        k <- Rinds2D[v,1]
		        l <- Rinds2D[v,2]
		        if (i==k) {
		            mH.mR.mR[u,v] <-  -mM[j,l]
				        if (i==j) {
				            mH.mR.mR[u,v] <- mH.mR.mR[u,v]*mR[i,i]
				        } 
				        if (k==l) {
				            mH.mR.mR[u,v] <- mH.mR.mR[u,v]*mR[k,k]
				        }		        
				        if ((i==j)&(j==k)&(k==l)) {
				            mH.mR.mR[u,v] <- mH.mR.mR[u,v] - mA[i,i]*mR[i,i]
				        }		            
		        }
        }
    }
    
    ###########################################################################
    
    mH.00 <- matrix(0,p+q,p+q)
    mH.00[1:p,1:p]         <-  -t(mX*vB2)%*%mX
    mH.00[p+(1:q),p+(1:q)] <-   mH.mR.mR
    
    ###########################################################################
    
    count <- 0
    linds <- list()
    for (i in 1:m) {
        linds[[i]] <- count + (1:nrow(lmZ[[i]]))
        count <- count + nrow(lmZ[[i]])
    }
    
    IK <- diag(1,K)
    
    lmH.0i <- list()
    lmH.ii <- list()
    for (i in 1:m) {
        inds <- linds[[i]]
        ni   <- length(inds)
        mXi  <- mX[inds,]
        mZi  <- lmZ[[i]]
        vmui <- vtheta[p+q+(i-1)*(K+q)+(1:K)]
        vri  <- vtheta[p+q+(i-1)*(K+q)+K+(1:q)]
        mRi[Rinds] <- vri
        mRi[Dinds] <- exp(mRi[Dinds])
        vB2i <- vB2[inds]
        vB3i <- vB3[inds]
        vB4i <- vB4[inds]
        
        #######################################################################
        
        mVi <- matrix(0,ni,q)
        for (j in 1:ni) {
            mA <- 2*mRi%*%mZi[j,]%*%t(mZi[j,]) 
            mA[Dinds] <- mA[Dinds]*mRi[Dinds]
            mVi[j,]   <- mA[Rinds]
        }
        
        #######################################################################
        
        vvi <- mR%*%vmui    
        mH.mR.vmui    <- matrix(0,q,K)
        for (s in 1:q) {
           j <- Rinds2D[s,1] 
           k <- Rinds2D[s,2] 
           if (j==k) {
               mH.mR.vmui[s,] <- -(mR[j,]*vmui[k] + IK[,j]*vvi[k])*mR[j,k]
           } else { 
               mH.mR.vmui[s,] <- -(mR[j,]*vmui[k] + IK[,k]*vvi[j])
           }
        }
        
        #######################################################################
        
        mA <- mRi%*%t(mR)
        mH.mR.mRi <- matrix(0,q,q)
		    for (u in 1:q) {
			      I <- Rinds2D[u,1]
				    J <- Rinds2D[u,2]		        
		        for (v in 1:q) {
				        j <- Rinds2D[v,1]
				        k <- Rinds2D[v,2]
				        mH.mR.mRi[u,v] <-  -mRi[j,J]*mR[I,k]
						    if (k==J) {
						        mH.mR.mRi[u,v] <- mH.mR.mRi[u,v] - mA[j,I]
						    } 
						    if (j==k) {
						        mH.mR.mRi[u,v] <- mH.mR.mRi[u,v]*mRi[j,j]
						    }
						    if (I==J) {
						        mH.mR.mRi[u,v] <- mH.mR.mRi[u,v]*mR[I,I]
						    }
				        #mH.mR.mR[v,u] <- mH.mR.mR[u,v]
		        }
		    }       
		    
		    #######################################################################
		     
        mH.0i <- matrix(0,p+q,K+q)
        mH.0i[1:p,1:K] <- -t(mXi*vB2i)%*%mZi
        mH.0i[1:p,K+(1:q)] <- -0.5*t(mXi*vB3i)%*%mVi
        mH.0i[p+(1:q),1:K] <- mH.mR.vmui
        mH.0i[p+(1:q),K+(1:q)] <- mH.mR.mRi
        lmH.0i[[i]] <- mH.0i
       
        #######################################################################
        
        mA0 <- -0.25*t(mVi*vB4i)%*%mVi
        mA1 <- (t(mZi*vB2i)%*%mZi + mSigma.inv) 
        mA2 <- mRi%*%mA1
        
        mH.ii.mRi.mRi <- matrix(0,q,q)
        for (u in 1:q) {
			      j <- Rinds2D[u,1]
				    k <- Rinds2D[u,2]		        
		        for (v in 1:q) {
				        J <- Rinds2D[v,1]
				        K <- Rinds2D[v,2]
				        mH.ii.mRi.mRi[u,v] <-   -as.numeric(j==J)*mA1[k,K] 								    
		        }
		    }
        for (u in 1:q) {
			      j <- Rinds2D[u,1]
				    k <- Rinds2D[u,2]		        
				    if (j==k) {
						    mH.ii.mRi.mRi[u,] <- mH.ii.mRi.mRi[u,]*mRi[j,j]
						    mH.ii.mRi.mRi[,u] <- mH.ii.mRi.mRi[,u]*mRi[j,j]
						    mH.ii.mRi.mRi[u,u] <- mH.ii.mRi.mRi[u,u] - mA2[j,j]*mRi[j,j]
				    }
		    }		    
		    mH.ii.mRi.mRi <- mH.ii.mRi.mRi + mA0
		    		    
		    #######################################################################
        
        mH.ii <- matrix(0,K+q,K+q)
        mH.ii[1:K,1:K] <-  -mA1
        mH.ii[1:K,K+(1:q)] <-  -0.5*t(mZi*vB3i)%*%mVi
        mH.ii[K+(1:q),1:K] <- t(mH.ii[1:K,K+(1:q)])
        mH.ii[K+(1:q),K+(1:q)] <- mH.ii.mRi.mRi
        lmH.ii[[i]] <- mH.ii
        
        #######################################################################
    }
    
    res.mH <- list(mH.00=mH.00,lmH.0i=lmH.0i,lmH.ii=lmH.ii)
            
		return(list(res.vg=res.vg,res.mH=res.mH))
}

###############################################################################

GVA.grouped.NRDIR <- function(vtheta,family,vy,mX,lmZ,gh)    
{
    m <- length(lmZ)
    
    res <- GVA.grouped.vgANDmH(vtheta,family,vy,mX,lmZ,gh)    
    res.vg <- res$res.vg
    res.mH <- res$res.mH
    vg <- c(res.vg$vg.0)
    vg.til <- res.vg$vg.0
    mH.til <- res.mH$mH.00
    for (i in 1:m) {
        vg.i <- res.vg$lvg.i[[i]]
        mH.0i     <- res.mH$lmH.0i[[i]]
        mH.ii.inv <- solve(res.mH$lmH.ii[[i]],tol=1.0E-99)
        mA <- mH.0i%*%mH.ii.inv
        vg <- c(vg,vg.i)
        vg.til <- vg.til - mA%*%vg.i
        mH.til <- mH.til - mA%*%t(mH.0i)
    } 

    sI <- solve(mH.til)
    vs <- sI%*%vg.til

    dvtheta <- vs
    for (i in 1:m) {
        dvtheta <- c(dvtheta,solve(res.mH$lmH.ii[[i]],res.vg$lvg.i[[i]] - t(res.mH$lmH.0i[[i]])%*%vs))
    }
  		
    return(list(dvtheta=-dvtheta,sI=-sI,vg=vg))
}

###############################################################################

getTheta <- function(vbeta,mSigma,lvmu,lmLambda) 
{    
    m <- length(lvmu)
    K <- ncol(mSigma)
    
    Dinds <- K*((1:K)-1) + (1:K) 	
    Rinds <- which(upper.tri(diag(0,K),diag=TRUE))	      
                      
    mSigma.inv <- solve(mSigma,tol=1.0E-99)  
    mR <- chol(mSigma.inv)
    mR[Dinds] <- log(mR[Dinds])

    vtheta <- c(vbeta,mR[Rinds])
    for (i in 1:m) {
        mR <- chol(lmLambda[[i]])    
		    mR[Dinds] <- log(mR[Dinds])
		    vtheta <- c(vtheta,lvmu[[i]],mR[Rinds])
    }
    return(vtheta)
}

###############################################################################

thetaToPars <- function(vtheta,vy,mX,lmZ) 
{
    n <- length(vy)
    m <- length(lmZ)
    K <- ncol(lmZ[[1]])
    
    Dinds <- K*((1:K)-1) + (1:K) 	
    Rinds <- which(upper.tri(diag(0,K),diag=TRUE))	  
    Linds <- which(lower.tri(diag(0,K),diag=TRUE))	  
    q <- length(Rinds)
        
    vbeta <- vtheta[1:p]
    vr    <- vtheta[p+(1:q)]
    mR <- matrix(0,K,K)
    mR[Rinds] <- vr
    mR[Dinds] <- exp(mR[Dinds])
    mSigma.inv <- t(mR)%*%mR   
    mSigma     <- solve(mSigma.inv,tol=1.0E-99) 
  
    mRi <- matrix(0,K,K)
    lvmu     <- list()
    lmLambda <- list()
    for (i in 1:m) {
        lvmu[[i]] <- vtheta[p+q+(i-1)*(K+q)+(1:K)]
        vri  <- vtheta[p+q+(i-1)*(K+q)+K+(1:q)]
        mRi[Rinds] <- vri
        mRi[Dinds] <- exp(mRi[Dinds])
        lmLambda[[i]] <- t(mRi)%*%mRi
    }
    
    return(list(vbeta=vbeta,mSigma=mSigma,lvmu=lvmu,lmLambda=lmLambda))
}

###############################################################################

GVA.grouped.FIT <- function(vbeta,mSigma,lvmu,lmLambda,family,vy,mX,lmZ) 
{
    MAXITER  <- 100
    EPS.TERM.PAR  <- 1.0E-5
    EPS.TERM.GRAD <- 1.0E-8
        
    gh <- NULL
    if (family=="LOGISTIC") {
				# Generate Quadrature Points (and cacluate w.til)
				N <- 15
				gh.temp <- gauss.quad(N,kind="hermite")
				gh <- list(x=gh.temp$nodes,w=gh.temp$weights,w.til=gh.temp$weights*exp(gh.temp$nodes^2))  
		}  
		
		vtheta <- getTheta(vbeta,mSigma,lvmu,lmLambda) 
    sl <- GVA.grouped.f(vtheta,family,vy,mX,lmZ,gh)  
    
    cat("Initial variational log-likelihood = ",sl)
    cat("\n\n")
    cat("ITER    VARLOGLIK        STEPSIZE         GRADERR       MAXPARDIF \n")
       
    for (ITER in 1:MAXITER) {
        vtheta.old <- vtheta
    
        # Calculate Newton-Raphson Direction 
        res <- GVA.grouped.NRDIR(vtheta,family,vy,mX,lmZ,gh)    
        
				stepsize <- 1                
        vtheta.temp <- vtheta + stepsize*res$dvtheta
        sl.temp <- GVA.grouped.f(vtheta.temp,family,vy,mX,lmZ,gh)  
        
        if (is.nan(sl.temp)) { 
            DOGRIDSEARCH <- TRUE
        } else {            
            if (sl.temp>sl) { DOGRIDSEARCH <- FALSE } 
            else { DOGRIDSEARCH <- TRUE }
        }         
        
        if (DOGRIDSEARCH) {
            NG <- 15
		        vsl <- c()
		        
		        # Allow for small negative step sizes 
		        # (just in case Hessian is negitive definite)
		        vstepsize <- seq(-0.25,1,,NG)[-NG]  
		        
		        for (JTER in 1:length(vstepsize)) {
		            vtheta.temp <- vtheta + vstepsize[JTER]*res$dvtheta
		            vsl[JTER]   <- GVA.grouped.f(vtheta.temp,family,vy,mX,lmZ,gh)  
		        }                   
		        sl <- max(vsl)
		        stepsize <- vstepsize[which.max(vsl)]
		        vtheta <- vtheta + stepsize*res$dvtheta
		    } else {
		        sl <- sl.temp
		        vtheta <- vtheta.temp    
        }
    
        # Check termination condtions
        err.par <- max(abs(vtheta-vtheta.old))/max(abs(vtheta))
        err.vg  <- sum(res$vg^2)/length(vtheta)
        
        cat(sprintf("%4d    %E    %+E    %E  %E\n",ITER,sl,stepsize,err.vg,err.par))

        if (err.vg<EPS.TERM.GRAD) { break; }
        if (err.par<EPS.TERM.PAR) { break; }
        
    }
    
    pars <- thetaToPars(vtheta,vy,mX,lmZ)
    vbeta.serr <- sqrt(diag(res$sI[1:ncol(mX),1:ncol(mX)]))
    
    return(list(vbeta=pars$vbeta,mSigma=pars$mSigma,lvmu=pars$lvmu,lmLambda=pars$lmLambda,vbeta.serr=vbeta.serr,sl=sl))
}

###############################################################################