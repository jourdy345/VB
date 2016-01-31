########## R script: Figure3EpilRISM.Rs #######################################
#
# Fit a poisson random slope model for the 
# Epilepsy Dataset using various methods.
#
# Last changed: 1 FEB 2011
#
###############################################################################

library(lme4)
library(nlme)
library(MASS)
load("epil2.rda")

###############################################################################

# Set various flags. 

do.PQL  <- TRUE
do.AGHQ <- TRUE
do.GVA  <- TRUE

###############################################################################

if (do.PQL) 
{
		cat("----------------------------------------------------------------------")
		cat("\n")
		
		cat("Fitting Random Slope Model for the Epilepsy Dataset \n")
		cat("using the PQL method from the MASS package\n") 
		cat("\n")
		
		bval1  <- proc.time()
    fit <- glmmPQL(y~Base*trt+Age+Visit,random=~Visit|subject,data=epil2,family="poisson")       
		eval1  <- proc.time()
		
		cat("Time taken = \n")
		cat(eval1[3]-bval1[3])
		cat("\n\n")		
				
		cat("PL estimates and confidence intervals are:\n")
		print(intervals(fit))
		
		cat("\n")
		cat("----------------------------------------------------------------------")	
		cat("\n")	   
}

###############################################################################

if (do.AGHQ) 
{
		cat("----------------------------------------------------------------------")
		cat("\n")
		
		cat("Fitting Random Slope Model for the Epilepsy Dataset \n")
		cat("using the AGHQ method from the lme4 package\n") 
		cat("\n")
		
		bval1  <- proc.time()
		aghq.fit <- glmer(formula=y~Base*trt+Age+Visit+(Visit|subject), data=epil2, family = poisson, nAGQ = 20)
		sum.aghq <- summary(aghq.fit)
		eval1  <- proc.time()
		
		cat("Time taken = \n")
		cat(eval1[3]-bval1[3])
		cat("\n\n")		
				
		print(sum.aghq)
		
		cat("\n")
		cat("----------------------------------------------------------------------")	
		cat("\n")			
}

###############################################################################

if (do.GVA) {		
		cat("----------------------------------------------------------------------")
		cat("\n")
		
		cat("Fitting Random Slope Model \n")
		cat("for the Epilepsy Dataset \n")
		cat("using the GVA method\n") 
		cat("\n")
    
		#############################################################################    
    
    nobs <- nrow(epil2)
				
    # set data matrices
		vy <- matrix(epil2$y,nobs,1)
		idNum <- epil2$subject
		mX <- c()
		lvy <- list()
	  lmX <- list()
		lmZ <- list()
		family <- "POISSON"
		
		subject <- 1
		vy.temp <- c()
		mX.temp <- c()
		mZ.temp <- c()		
		m <- 1
		for (i in 1:nobs) {
				if (epil2$subject[i]!=subject) {
						lvy[[m]] <- vy.temp
				    lmX[[m]] <- mX.temp 	
				    lmZ[[m]] <- mZ.temp
					  rownames(lmX[[m]]) <- NULL
					  rownames(lmZ[[m]]) <- NULL				    
						vy.temp <- c()
						mX.temp <- c()
						mZ.temp <- c()		
						m <- m + 1
						subject <- subject + 1
		    }
		    trt <- as.numeric(epil2$trt[i]!="placebo")
				vx <- c(1,epil2$Base[i],trt,epil2$Base[i]*trt,epil2$Age[i],epil2$Visit[i])
				vz <- c(1,epil2$Visit[i])
				mX <- rbind(mX,vx)		
				vy.temp <- rbind(vy.temp,epil2$y[i])		
				mX.temp <- rbind(mX.temp,vx)
				mZ.temp <- rbind(mZ.temp,vz)		    
		}		
		lvy[[m]] <- vy.temp
		lmX[[m]] <- mX.temp 	
	  lmZ[[m]] <- mZ.temp	
	  rownames(lmX[[m]]) <- NULL
	  rownames(lmZ[[m]]) <- NULL	  
	  rownames(mX) <- NULL	

		m    <- 59
		n    <- 4
		p    <- ncol(lmX[[1]])
		K    <- ncol(lmZ[[1]])    
		
		#############################################################################
    
    vbeta  <- matrix(0,p,1)
    mSigma <- diag(1,K)
    lvmu <- list()
    lmLambda <- list()
    for (i in 1:m) {
        lvmu[[i]]     <- matrix(0,K,1)
        lmLambda[[i]] <- diag(1,K)
    }
  				    
		bval1  <- proc.time()				
		source("GVAgrouped.Rs")				
		resVar  <- GVA.grouped.FIT(vbeta,mSigma,lvmu,lmLambda,family,vy,mX,lmZ) 
		eval1  <- proc.time()
		print(eval1[3]-bval1[3])
		
		for (i in 1:p) {
		    cat(sprintf("beta[%1d] = %+E (%E) \n",i-1,resVar$pars$vbeta[i],resVar$vbeta.serr[i]))
		}
		
		sigma0 <- sqrt(resVar$pars$mSigma[1,1])
		sigma1 <- sqrt(resVar$pars$mSigma[2,2])
		rho <- resVar$pars$mSigma[1,2]/(sigma0*sigma1)
		
		cat("sigma0 = ",sigma0,"\n")
		cat("sigma1 = ",sigma1,"\n")
		cat("rho = ",rho,"\n")
		
		
		cat("\n")
		cat("----------------------------------------------------------------------")	
		cat("\n")			
}
		  
########## End of Figure3EpilRISM.Rs ##########################################
