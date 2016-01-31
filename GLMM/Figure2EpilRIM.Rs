########## R script: Figure2EpilRIM.Rs ##########

# Fit a poisson random intercept model 
# for the Epilepsy Dataset
# using various methods.

# Last changed: 30 NOV 2010

###############################################################################

library(lme4)
library(MASS)

source("GVArandint.Rs")
load("epil2.rda")

###############################################################################

# Set various flags. 
# 
do.PQL  <- TRUE
do.AGHQ <- TRUE
do.GVA  <- TRUE

###############################################################################

if (do.PQL) 
{
		cat("----------------------------------------------------------------------")
		cat("\n")
		
		cat("Fitting Random Intercept Model \n")
		cat("for the Epilepsy Dataset \n")
		cat("using the PQL method from the MASS package\n") 
		cat("\n")
		
		bval1  <- proc.time()
    fit <- glmmPQL(y~Base*trt+Age+V4,random=~1|subject,data=epil2,family=poisson)
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
		cat("for the Epilepsy Dataset \n")
		cat("using the AGHQ method from the lme4 package\n") 
		cat("\n")
		
		bval1  <- proc.time()
		aghq.fit <- glmer(formula=y~Base*trt+Age+V4+(1|subject), data=epil2, family = poisson, nAGQ = 20)
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

    # Set up problem for GVA to solve
		vy    <- epil2$y
		n     <- length(vy)
		trt  <- as.numeric(epil2$trt!="placebo")
		mX   <- cbind(1,epil2$Base,trt,epil2$Age,epil2$V4,epil2$Base*trt)
		p    <- ncol(mX)
		id   <- epil2$subject
		m    <- length(unique(id))	
		vn   <- rep(4,m)
		
		# Initialise variables
		vbeta <- matrix(0,p,1)
		gamma <- -1
		vmu   <- matrix(0,m,1)
		vzeta <- matrix(-2,m,1)		
				
		#############################################################################
				
		# Perform Fitting
		bval1  <- proc.time()
		res <- GVA.randint.FIT(vbeta,exp(gamma),vmu,exp(vzeta),family="POISSON",vy,mX,vn,id) 
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

