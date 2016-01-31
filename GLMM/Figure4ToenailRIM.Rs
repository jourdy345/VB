########## R script: Figure2ToenailRIM.Rs ##########

# Fit a poisson random intercept model 
# for the Toenail Dataset
# using various methods.

# Last changed: 30 NOV 2010

###############################################################################

rm(list = ls())
memory.limit(size=2048)

source("GVArandint.Rs")

library(lme4)
library(MASS)
library(nlme)
load("toenail.rda")

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
		
		cat("Fitting Random Intercept Model \n")
		cat("for the Toenail Dataset \n")
		cat("using the PQL method from the MASS package\n") 
		cat("\n")
		
		bval1  <- proc.time()
    fit <- glmmPQL(infect~trt*time,data=toenail,random=~1|idnr,family=binomial) 
		sum.fit <- summary(fit)
		sigma.pql.hat <- exp(as.numeric(fit$modelStruct$reStruct))			        
		eval1  <- proc.time()
		
		cat("Time taken = \n")
		cat(eval1[3]-bval1[3])
		cat("\n\n")		
				
		cat("Estimated Coefficients = \n")
		cat( sprintf("%+.3f",round(fit$coefficients$fixed,3)), "\n" )
		cat("\n")
		
		cat("Estimated Coefficient Standard Errors = \n")
		cat( sprintf("%+.3f",round(sum.fit$tTable[,2],3)), "\n" )
		cat("\n")
		
		cat("Estimated Variance Component Value = \n")
		cat(round(sigma.pql.hat^2,3))
		
		cat("\n")
		cat("----------------------------------------------------------------------")	
		cat("\n")	    
}

###############################################################################

if (do.AGHQ) 
{
		cat("----------------------------------------------------------------------")
		cat("\n")
		
		cat("Fitting Random Intercept Model \n")
		cat("for the Toenail Dataset \n")
		cat("using the AGHQ method from the lme4 package\n") 
		cat("\n")
		
		bval1  <- proc.time()
		aghq.fit <- glmer(formula=infect~trt*time + (1|idnr), data=toenail, family = binomial, nAGQ = 30)
		sum.aghq <- summary(aghq.fit)
		eval1  <- proc.time()
		
		cat("Time taken = \n")
		cat(eval1[3]-bval1[3])
		cat("\n\n")		
				
		cat("Estimated Coefficients = \n")
		cat( sprintf("%+.3f",round(aghq.fit@fixef,3)), "\n" )
		cat("\n")
		
		cat("Estimated Coefficient Standard Errors = \n")
		vbeta.serr <- attributes(sum.aghq)$coef[,2]
		cat( sprintf("%+.3f",round(vbeta.serr,3)), "\n" )
		cat("\n")
		
		cat("Estimated Variance Component Value = \n")
		cat(round(as.numeric(sum.aghq@REmat[3]),3))
		
		cat("\n")
		cat("----------------------------------------------------------------------")	
		cat("\n")			
}

###############################################################################

if (do.GVA) 
{
		cat("----------------------------------------------------------------------")
		cat("\n")
		
		cat("Fitting Random Intercept Model \n")
		cat("for the Epilepsy Dataset \n")
		cat("using the GVA method\n") 
		cat("\n")
		
		vy <- toenail$infect
		n <- length(vy)
		
		mX <- cbind(1,toenail$trt,toenail$time,toenail$time*toenail$trt)
    p <- ncol(mX)
    
    unique.id <- unique(toenail$idnr)
    m <- length(unique.id)
    
    id <- c()
    vn <- c()
    for (i in 1:m) { 
        ind  <- which(toenail$idnr==unique.id[i])
        vn[i] <- length(ind)
        id[ind] <- i
    }
    		
		###############################################################################
		
		# Initialise variables
		vbeta <- matrix(0,p,1) #matrix(aghq.fit@fixef)
		gamma <- 2             #as.numeric(sum.aghq@REmat[4])
		vmu   <- matrix(0,m,1) #matrix(aghq.fit@ranef)
		vzeta <- matrix(2,m,1)
		
		###############################################################################
		
		# Perform Fitting
		bval1  <- proc.time()
		res <- GVA.randint.FIT(vbeta,exp(gamma),vmu,exp(vzeta),family="LOGISTIC",vy,mX,vn,id) 
		eval1  <- proc.time()

		#######################################
		
		cat("Time taken = \n")
		cat(eval1[3]-bval1[3])
		cat("\n\n")
		
		cat("Variational log-likeihood = \n")
		cat(res$sl)
		cat("\n\n")
		
		cat("Estimated Coefficients = \n")
		cat( sprintf("%+.3f",round(res$vbeta,3)), "\n" )
		cat("\n")
		
		cat("Estimated Coefficient Standard Errors = \n")
		cat( sprintf("%+.3f",round(res$vbeta.serr,3)), "\n" )
		cat("\n")
		
		cat("Estimated Variance Component Value = \n")
		cat(round(res$sigma2,3))
		
		cat("\n")
		cat("----------------------------------------------------------------------")
		cat("\n")
}

###############################################################################
