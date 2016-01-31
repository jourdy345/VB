########## R script: Figure2EpilRIM.Rs ##########

# Fit a poisson random intercept model 
# for the Epilepsy Dataset
# using various methods.

# Last changed: 30 NOV 2010

###############################################################################

#rm(list = ls())
memory.limit(size=2048)

source("GVArandint.Rs")

library(lme4)
library(MASS)

do.PQL  <- TRUE
do.AGHQ <- TRUE
do.LAP  <- TRUE
do.GVA  <- TRUE

###############################################################################

trapint <- function(xgrid, fgrid) 
{
    ng <- length(xgrid)
    xvec <- xgrid[2:ng] - xgrid[1:(ng - 1)]
    fvec <- fgrid[1:(ng - 1)] + fgrid[2:ng]
    integ <- sum(xvec * fvec)/2
    return(integ)
}

###############################################################################

randEffDist <- function(vy,mX,vbeta,sigma2,mu,lambda) 
{
    m <- length(id) 
    u.g <- seq(mu-5*sqrt(lambda),mu+5*sqrt(lambda),,1000)
    sl.g <- 0*u.g
    for (i in 1:length(u.g)) {
        veta  <- mX%*%vbeta + u.g[i]
        sl.g[i] <- sum(vy*veta - log(1+exp(veta))) - 0.5*u.g[i]^2/sigma2
    }
    sl.g <- exp(sl.g - max(sl.g))
    vf.g <- sl.g/trapint(u.g,sl.g) 
    mu     <- trapint(u.g,u.g*vf.g)
    lambda <- trapint(u.g,(u.g - mu)^2*vf.g)
    skew   <- trapint(u.g,(u.g - mu)^3*vf.g)/(lambda^(3/2))
    return(list(x=u.g,y=vf.g,mu=mu,lambda=lambda,skew=skew))
}

###############################################################################

lap.randEffDist <- function(vy,mX,id,ni,vbeta,sigma2,vmu,vlambda) 
{
    MAXITER <- 100
    m <- length(id)    
    for (ITER in 1:MAXITER) {                     
        veta <- mX%*%vbeta + rep(vmu,each=ni)
        mu   <- 1/(1+exp(-veta))
        veps <- vy - mu
        vw   <- mu*(1 - mu)
        sw   <- tapply(vw,id,sum)
        seps <- tapply(veps,id,sum)
        vlambda <- 1/(sw + 1/sigma2)
        dvmu <- vlambda*( seps - vmu/sigma2 )
        vmu  <- vmu + 0.5*dvmu
        verr <- abs(dvmu)
        if (max(verr) < 5.0e-5) { break; }
    }  
    return(list(vmu=vmu,vlambda=vlambda)) 
}

###############################################################################

gva.randEffDist <- function(vy,mX,id,ni,vbeta,sigma2,vmu,vlambda,gh) 
{
		# Generate Quadrature Points
		N <- 15
		gh.temp <- gauss.quad(N,kind="hermite")
		gh <- list(x=gh.temp$nodes,w=gh.temp$weights,w.til=gh.temp$weights*exp(gh.temp$nodes^2))  
 
    MAXITER <- 100
    m <- length(id)    
    for (ITER in 1:MAXITER) {                     
        veta.vmu     <- mX%*%vbeta + rep(vmu,each=ni)
        veta.vsigma2 <- rep(vlambda,each=ni)
        res <- B12.fun("LOGISTIC",veta.vmu,veta.vsigma2,gh) 
        sw   <- tapply(res$vB2,id,sum)
        seps <- tapply(vy - res$vB1,id,sum)
        vlambda <- 1/(sw + 1/sigma2)
        dvmu <- vlambda*( seps - vmu/sigma2 )
        vmu  <- vmu + 0.5*dvmu
        verr <- abs(dvmu)
        if (max(verr) < 5.0e-5) { break; }
    }  
    return(list(vmu=vmu,vlambda=vlambda)) 
}

###############################################################################

# Settings
SIGMA <- sqrt(1.5)
M     <- c(15,50)
ni    <- 8
beta0.true <- 0
beta1.true <- 5

TRIALS  <- 2000  # Total number of trials
TRIALS2 <- 200   # Trials where u|y is calculated exactly

###############################################################################

# Set up tables to record results

beta0.pql  <- matrix(0,length(M)*length(SIGMA)*TRIALS,1); dim(beta0.pql)  <- c(length(M),length(SIGMA),TRIALS)
beta0.lap  <- matrix(0,length(M)*length(SIGMA)*TRIALS,1); dim(beta0.lap)  <- c(length(M),length(SIGMA),TRIALS)
beta0.aghq <- matrix(0,length(M)*length(SIGMA)*TRIALS,1); dim(beta0.aghq) <- c(length(M),length(SIGMA),TRIALS)
beta0.gva  <- matrix(0,length(M)*length(SIGMA)*TRIALS,1); dim(beta0.gva)  <- c(length(M),length(SIGMA),TRIALS)

beta1.pql  <- matrix(0,length(M)*length(SIGMA)*TRIALS,1); dim(beta1.pql)  <- c(length(M),length(SIGMA),TRIALS)
beta1.lap  <- matrix(0,length(M)*length(SIGMA)*TRIALS,1); dim(beta1.lap)  <- c(length(M),length(SIGMA),TRIALS)
beta1.aghq <- matrix(0,length(M)*length(SIGMA)*TRIALS,1); dim(beta1.aghq) <- c(length(M),length(SIGMA),TRIALS)
beta1.gva  <- matrix(0,length(M)*length(SIGMA)*TRIALS,1); dim(beta1.gva)  <- c(length(M),length(SIGMA),TRIALS)

sigma.pql  <- matrix(0,length(M)*length(SIGMA)*TRIALS,1); dim(sigma.pql)  <- c(length(M),length(SIGMA),TRIALS)
sigma.lap  <- matrix(0,length(M)*length(SIGMA)*TRIALS,1); dim(sigma.lap)  <- c(length(M),length(SIGMA),TRIALS)
sigma.aghq <- matrix(0,length(M)*length(SIGMA)*TRIALS,1); dim(sigma.aghq) <- c(length(M),length(SIGMA),TRIALS)
sigma.gva  <- matrix(0,length(M)*length(SIGMA)*TRIALS,1); dim(sigma.gva)  <- c(length(M),length(SIGMA),TRIALS)

stderr.beta0.pql  <- matrix(0,length(M)*length(SIGMA)*TRIALS,1); dim(stderr.beta0.pql) <- c(length(M),length(SIGMA),TRIALS)
stderr.beta0.lap  <- matrix(0,length(M)*length(SIGMA)*TRIALS,1); dim(stderr.beta0.lap) <- c(length(M),length(SIGMA),TRIALS)
stderr.beta0.aghq <- matrix(0,length(M)*length(SIGMA)*TRIALS,1); dim(stderr.beta0.aghq) <- c(length(M),length(SIGMA),TRIALS)
stderr.beta0.gva  <- matrix(0,length(M)*length(SIGMA)*TRIALS,1); dim(stderr.beta0.gva) <- c(length(M),length(SIGMA),TRIALS)

stderr.beta1.pql  <- matrix(0,length(M)*length(SIGMA)*TRIALS,1); dim(stderr.beta1.pql) <- c(length(M),length(SIGMA),TRIALS)
stderr.beta1.lap  <- matrix(0,length(M)*length(SIGMA)*TRIALS,1); dim(stderr.beta1.lap) <- c(length(M),length(SIGMA),TRIALS)
stderr.beta1.aghq <- matrix(0,length(M)*length(SIGMA)*TRIALS,1); dim(stderr.beta1.aghq) <- c(length(M),length(SIGMA),TRIALS)
stderr.beta1.gva  <- matrix(0,length(M)*length(SIGMA)*TRIALS,1); dim(stderr.beta1.gva) <- c(length(M),length(SIGMA),TRIALS)

times.pql  <- matrix(0,length(M)*length(SIGMA)*TRIALS,1); dim(times.pql) <- c(length(M),length(SIGMA),TRIALS)
times.lap  <- matrix(0,length(M)*length(SIGMA)*TRIALS,1); dim(times.lap) <- c(length(M),length(SIGMA),TRIALS)
times.aghq <- matrix(0,length(M)*length(SIGMA)*TRIALS,1); dim(times.aghq) <- c(length(M),length(SIGMA),TRIALS)
times.gva  <- matrix(0,length(M)*length(SIGMA)*TRIALS,1); dim(times.gva) <- c(length(M),length(SIGMA),TRIALS)

converged.pql  <- matrix(0,length(M)*length(SIGMA)*TRIALS,1); dim(converged.pql) <- c(length(M),length(SIGMA),TRIALS)
converged.lap  <- matrix(0,length(M)*length(SIGMA)*TRIALS,1); dim(converged.lap) <- c(length(M),length(SIGMA),TRIALS)
converged.aghq <- matrix(0,length(M)*length(SIGMA)*TRIALS,1); dim(converged.aghq) <- c(length(M),length(SIGMA),TRIALS)
converged.gva  <- matrix(0,length(M)*length(SIGMA)*TRIALS,1); dim(converged.gva) <- c(length(M),length(SIGMA),TRIALS)

rmse.vmu.lap  <- matrix(0,length(M)*length(SIGMA)*TRIALS2,1); dim(rmse.vmu.lap) <- c(length(M),length(SIGMA),TRIALS2)
rmse.vmu.gva  <- matrix(0,length(M)*length(SIGMA)*TRIALS2,1); dim(rmse.vmu.gva) <- c(length(M),length(SIGMA),TRIALS2)

###############################################################################

ITER0 <- 1
JTER0 <- 1
trial0 <-  1  

for (ITER in ITER0:length(M)) {
    ITER0  <- 1				
    for (JTER in JTER0:length(SIGMA)) {
				JTER0  <- 1				
				for (trial in trial0:TRIALS) {
				    trial0 <- 1        
				
				    #############################################################################
				    
				    # Simulate Data

				    m <- M[ITER]
				    sigma.true <- SIGMA[JTER]
				    
				    set.seed(trial)
				    
				    vu.true <- rnorm(m,0,sigma.true)
				    veta <- c()
				    vx   <- c()
				    id   <- c()
				    for (i in 1:m) {
				        for (j in 1:ni) {
				            id <- c(id,i)
				            x <- j/ni
				            vx <- c(vx,x)
				            veta <- c(veta,beta0.true + beta1.true*x + vu.true[i])
				        }
				    }
				    vmu.true <- 1/(1+exp(-veta))
				    vy <- matrix(rbinom(ni*m,1,vmu.true))
				    				    
				    #############################################################################
				
						if (do.PQL) {
						    cat("trial=",trial,"    m=",m,"    sigma2=",sigma.true^2,"\n")
						    
						    cat("Fitting model using PQL \n")
						    
								bval1  <- proc.time()
								res <- try(fit <- glmmPQL(vy~vx,random=~1|id,family=binomial),silent=TRUE)
								eval1  <- proc.time()
								
								if ((class(res)[1])!="try-error") {
										sum.fit <- summary(fit)
										sigma.pql.hat <- exp(as.numeric(fit$modelStruct$reStruct))
										    
										print(c(fit$coefficients$fixed,sigma.pql.hat^2))
										print(eval1[3]-bval1[3])  
										  										    
										beta0.pql[ITER,JTER,trial] <- fit$coefficients$fixed[1]
										beta1.pql[ITER,JTER,trial] <- fit$coefficients$fixed[2]
										sigma.pql[ITER,JTER,trial] <- sigma.pql.hat
										stderr.beta0.pql[ITER,JTER,trial] <- sum.fit$tTable[1,2]
										stderr.beta1.pql[ITER,JTER,trial] <- sum.fit$tTable[2,2]
				            times.pql[ITER,JTER,trial] <- eval1[3]-bval1[3]
				            converged.pql[ITER,JTER,trial] <- TRUE
		            } else {
		                converged.pql[ITER,JTER,trial] <- FALSE
		            }
						}
						
						#############################################################################
						
						if (do.LAP) {
						    
						    last.warning <- NULL
						    
						    cat("trial=",trial,"    m=",m,"    sigma2=",sigma.true^2,"\n")
						    print("Laplace")
								bval1  <- proc.time()
								res <- try(lap.fit <- glmer(formula=vy~vx + (1|id), family = binomial, nAGQ = 1),silent=TRUE)
								sum.lap <- summary(lap.fit)
								eval1  <- proc.time()
								
								if ((class(res)[1])!="try-error") {
										print(c(lap.fit@fixef,as.numeric(sum.lap@REmat[3])))
										print(eval1[3]-bval1[3])
										
										beta0.lap[ITER,JTER,trial] <- lap.fit@fixef[1]
										beta1.lap[ITER,JTER,trial] <- lap.fit@fixef[2]
										sigma.lap[ITER,JTER,trial] <- as.numeric(sum.lap@REmat[4])
										stderr.beta0.lap[ITER,JTER,trial] <- sum.lap@coefs[1,2]
										stderr.beta1.lap[ITER,JTER,trial] <- sum.lap@coefs[2,2]								
										times.lap[ITER,JTER,trial] <- eval1[3]-bval1[3]		
										
										# Note that glmer sometimes gives a warning
										# suggesting convergence problems. If this 
										# is true we should doubt that glmer actually
										# converged.
				            
				            converged.lap[ITER,JTER,trial] <- is.null(last.warning)
		            } else {
		                converged.lap[ITER,JTER,trial] <- FALSE
		            }		            
						}
						
						#############################################################################
						
					  if (do.AGHQ) {
					  
					      last.warning <- NULL
					  
					      cat("trial=",trial,"    m=",m,"    sigma2=",sigma.true^2,"\n")
					  
						    print("AGHQ")
								bval1  <- proc.time()
								res <- try(aghq.fit <- glmer(formula=vy~vx + (1|id), family = binomial, nAGQ = 15),silent=TRUE)
								sum.aghq <- summary(aghq.fit)
								eval1  <- proc.time()		
								
								if ((class(res)[1])!="try-error") {								
										print(c(aghq.fit@fixef,as.numeric(sum.aghq@REmat[3])))		
										print(eval1[3]-bval1[3])
										
										beta0.aghq[ITER,JTER,trial] <- aghq.fit@fixef[1]
										beta1.aghq[ITER,JTER,trial] <- aghq.fit@fixef[2]
										sigma.aghq[ITER,JTER,trial] <- as.numeric(sum.aghq@REmat[4])
										stderr.beta0.aghq[ITER,JTER,trial] <- sum.aghq@coefs[1,2]
										stderr.beta1.aghq[ITER,JTER,trial] <- sum.aghq@coefs[2,2]
		
				            times.aghq[ITER,JTER,trial] <- eval1[3]-bval1[3]	
				            
										# Note that glmer sometimes gives a warning
										# suggesting convergence problems. If this 
										# is true we should doubt that glmer actually
										# converged.
				            
				            converged.aghq[ITER,JTER,trial] <- is.null(last.warning)     
		            } else {
		                converged.aghq[ITER,JTER,trial] <- FALSE
		            }			                   									
						}
						
						#############################################################################
						
						# Put data into form used by GVA
								
						mX <- cbind(1,vx)																		
						vn <- rep(ni,m)		

						###############################################################################
						
						# Set initial values
						
						vbeta <- matrix(0,p,1)  
						gamma <- log(2*sigma.true^2)  
						vmu   <- matrix(0,m,1)  
					
						
						
						if ((ITER==1)&(JTER==1)) { 
						    # Initial values for vzeta (which work most of the time)
						    vzeta <- matrix(2,m,1)
						    
						    # Some exceptions where GVA had problems so use different initial values
						   	indExep <- c(52,109,185,298,656,1074,1283,1296,1360,1574,1813,1845)
						    if (any(trial==indExep)) {
						        vzeta <- matrix( 1,m,1)
						    }     		    
						   	indExep <- c(1180)
						    if (any(trial==indExep)) {
						        vzeta <- matrix( 0.5,m,1)
						    }     						    					    			      					      
						}
						
						if ((ITER==2)&(JTER==1)) { 
                # Initial values for vzeta (which work most of the time)
						    vzeta <- matrix( 1,m,1)
						    
						    # Some exceptions where GVA had problems so use different initial values
						   	indExep <- c(1292)
						    if (any(trial==indExep)) {
						        vzeta <- matrix(0,m,1)
						    }						    				    					    					    					    						    					    							    					    
						}
						
						###############################################################################
						
						cat("trial=",trial,"    m=",m,"    sigma2=",sigma.true^2,"\n")
						
						print("GVA.randint.NR") 
						bval1  <- proc.time()
						res <- GVA.randint.FIT(vbeta,exp(gamma),vmu,exp(vzeta),family="LOGISTIC",vy,mX,vn,id) 
						eval1  <- proc.time()
						
						print(res$sl)
						print(round(c(res$vbeta,res$sigma2),3))
						print(eval1[3]-bval1[3])
						
						beta0.gva[ITER,JTER,trial] <- res$vbeta[1]
						beta1.gva[ITER,JTER,trial] <- res$vbeta[2]
						sigma.gva[ITER,JTER,trial] <- sqrt(res$sigma2)
						stderr.beta0.gva[ITER,JTER,trial] <- res$vbeta.serr[1]
						stderr.beta1.gva[ITER,JTER,trial] <- res$vbeta.serr[2]				
		        times.gva[ITER,JTER,trial] <- eval1[3]-bval1[3]	
		        converged.gva[ITER,JTER,trial] <- TRUE
						
						###############################################################################
						
						# Only calculate E(u|y) and Var(u|y) for the first 200 trials
						if (trial<=TRIALS2) {
						
						    # Setup matrices
						
								lvy <- list()    
								lmX <- list()
								ids <- unique(id) 
								for (i in 1:m) {		
								    ind <- (id==ids[i])
								    lmX[[i]] <- cbind(rep(1,sum(ind)),vx[ind])
								    lvy[[i]] <- vy[ind]
								}						
						
								vbeta.true  <- matrix(c(beta0.true,beta1.true))
								sigma2.true <- sigma.true^2
								vmu.init <- res$vmu
								vlambda.init <- res$vlambda
								
								# Using the true values calculate u|y using LAP and GVA methods
								res.lap <- lap.randEffDist(vy,mX,id,ni,vbeta.true,sigma2.true,vmu.init,vlambda.init) 
		            res.gva <- gva.randEffDist(vy,mX,id,ni,vbeta.true,sigma2.true,vmu.init,vlambda.init) 
		            
		            # Using the true values calculate u|y using numerical quadrature
		            vmu.true <- c()
		            for (i in 1:m) {
				            res.true <- randEffDist(matrix(lvy[[i]]),lmX[[i]],vbeta.true,sigma2.true,res.gva$vmu[i],res.gva$vlambda[i]) 
				            vmu.true[i]     <-  res.true$mu
				        }
								rmse.vmu.lap[ITER,JTER,trial]  <- sqrt( mean((res.lap$vmu - vmu.true)^2)	)
								rmse.vmu.gva[ITER,JTER,trial]  <- sqrt( mean((res.gva$vmu - vmu.true)^2)	)
						}
						
						###############################################################################					
						
						save.image("Table1Setting3.Rdata")	
				
				}
		}
}

###############################################################################

indOutliers <- function(x,nsd=5) {
    z.scores <- (x - median(x))/mad(x)
    return(which(abs(z.scores)>nsd))
}


for (ITER in 1:length(M)) {
     for (JTER in 1:length(SIGMA)) {
          
 				 m <- M[ITER]
 				 sigma.true <- SIGMA[JTER]    
 
         cat("m=",m,"    sigma=",sigma.true,"\n")
         
         # Remove cases where any one of the methods did not converge         
         is.included <- which( (converged.pql[ITER,JTER,]&converged.lap[ITER,JTER,]&converged.aghq[ITER,JTER,]&converged.gva[ITER,JTER,]))
         inds <- is.included
         
         # When sigma estimate is an outlier (probably due to algorith 
         # divergence) remove it from cases used in the analysis. 
         # We assume this is a property of the software rather than 
         # a property of the estimators themselves.
         # This usually makes LAP and AGHQ methods appear better.
         
         inds.outliers <- c(indOutliers(beta0.pql[ITER,JTER,inds]),
                            indOutliers(beta0.lap[ITER,JTER,inds]),
                            indOutliers(beta0.aghq[ITER,JTER,inds]),
                            indOutliers(beta0.gva[ITER,JTER,inds]),
														indOutliers(beta1.pql[ITER,JTER,inds]),
                            indOutliers(beta1.lap[ITER,JTER,inds]),
                            indOutliers(beta1.aghq[ITER,JTER,inds]),
                            indOutliers(beta1.gva[ITER,JTER,inds]),    
														indOutliers(sigma.pql[ITER,JTER,inds]),
                            indOutliers(sigma.lap[ITER,JTER,inds]),
                            indOutliers(sigma.aghq[ITER,JTER,inds]),
                            indOutliers(sigma.gva[ITER,JTER,inds]))                                                   
         if (length(inds.outliers)>0) {              
             inds <- inds[-inds.outliers]     
         }
         
         print(length(inds))

         locfun <- mean
         sprfun <- sd
         
         ######################################################################
         
         # Calculate and print summary statistic for fitted beta0 values
         
         cat("beta0 \n")
         
         dp <- 2
         
         mean.beta0.pql  <- round(locfun(beta0.pql[ITER,JTER,inds]),dp)
         mean.beta0.lap  <- round(locfun(beta0.lap[ITER,JTER,inds]),dp)
         mean.beta0.aghq <- round(locfun(beta0.aghq[ITER,JTER,inds]),dp)
         mean.beta0.gva  <- round(locfun(beta0.gva[ITER,JTER,inds]),dp)
         
         estderr.beta0.pql  <- round(sprfun(beta0.pql[ITER,JTER,inds])/sqrt(length(inds)),dp)
         estderr.beta0.lap  <- round(sprfun(beta0.lap[ITER,JTER,inds])/sqrt(length(inds)),dp)
         estderr.beta0.aghq <- round(sprfun(beta0.aghq[ITER,JTER,inds]/sqrt(length(inds))),dp)
         estderr.beta0.gva  <- round(sprfun(beta0.gva[ITER,JTER,inds]/sqrt(length(inds))),dp)           
         
         sd.beta0.pql  <- round(sprfun(beta0.pql[ITER,JTER,inds]),dp)
         sd.beta0.lap  <- round(sprfun(beta0.lap[ITER,JTER,inds]),dp)
         sd.beta0.aghq <- round(sprfun(beta0.aghq[ITER,JTER,inds]),dp)
         sd.beta0.gva  <- round(sprfun(beta0.gva[ITER,JTER,inds]),dp)  
         
         mstderr.beta0.pql  <- round(locfun(stderr.beta0.pql[ITER,JTER,inds]),dp)
         mstderr.beta0.lap  <- round(locfun(stderr.beta0.lap[ITER,JTER,inds]),dp)
         mstderr.beta0.aghq <- round(locfun(stderr.beta0.aghq[ITER,JTER,inds]),dp)
         mstderr.beta0.gva  <- round(locfun(stderr.beta0.gva[ITER,JTER,inds]),dp)
         
         rmse.beta0.pql  <- round(sqrt(locfun((beta0.pql[ITER,JTER,inds] - beta0.true)^2)),dp)
         rmse.beta0.lap  <- round(sqrt(locfun((beta0.lap[ITER,JTER,inds] - beta0.true)^2)),dp)
         rmse.beta0.aghq <- round(sqrt(locfun((beta0.aghq[ITER,JTER,inds] - beta0.true)^2)),dp)
         rmse.beta0.gva  <- round(sqrt(locfun((beta0.gva[ITER,JTER,inds] - beta0.true)^2)),dp)   
                  
         cat(mean.beta0.pql, "(",estderr.beta0.pql ,") & " ,"(",sd.beta0.pql ,") & ", mstderr.beta0.pql ," & ", rmse.beta0.pql, " \n")
         cat(mean.beta0.lap, "(",estderr.beta0.lap ,") & " ,"(",sd.beta0.lap ,") & ", mstderr.beta0.lap ," & ", rmse.beta0.lap, " \n")
         cat(mean.beta0.aghq,"(",estderr.beta0.aghq ,") & " ,"(",sd.beta0.aghq,") & ", mstderr.beta0.aghq," & ", rmse.beta0.aghq, " \n")
         cat(mean.beta0.gva, "(",estderr.beta0.gva ,") & " ,"(",sd.beta0.gva ,") & ", mstderr.beta0.gva ," & ", rmse.beta0.gva, " \n\n")
         
         ######################################################################
                  
         # Calculate and print summary statistic for fitted beta1 values
         
         cat("beta1 \n")
         
         mean.beta1.pql  <- round(locfun(beta1.pql[ITER,JTER,inds]),dp)
         mean.beta1.lap  <- round(locfun(beta1.lap[ITER,JTER,inds]),dp)
         mean.beta1.aghq <- round(locfun(beta1.aghq[ITER,JTER,inds]),dp)
         mean.beta1.gva  <- round(locfun(beta1.gva[ITER,JTER,inds]),dp)
         
         estderr.beta1.pql  <- round(sprfun(beta1.pql[ITER,JTER,inds])/sqrt(length(inds)),dp)
         estderr.beta1.lap  <- round(sprfun(beta1.lap[ITER,JTER,inds])/sqrt(length(inds)),dp)
         estderr.beta1.aghq <- round(sprfun(beta1.aghq[ITER,JTER,inds]/sqrt(length(inds))),dp)
         estderr.beta1.gva  <- round(sprfun(beta1.gva[ITER,JTER,inds]/sqrt(length(inds))),dp)             
         
         sd.beta1.pql  <- round(sprfun(beta1.pql[ITER,JTER,inds]),dp)
         sd.beta1.lap  <- round(sprfun(beta1.lap[ITER,JTER,inds]),dp)
         sd.beta1.aghq <- round(sprfun(beta1.aghq[ITER,JTER,inds]),dp)
         sd.beta1.gva  <- round(sprfun(beta1.gva[ITER,JTER,inds]),dp)  
         
         mstderr.beta1.pql  <- round(locfun(stderr.beta1.pql[ITER,JTER,inds]),dp)
         mstderr.beta1.lap  <- round(locfun(stderr.beta1.lap[ITER,JTER,inds]),dp)
         mstderr.beta1.aghq <- round(locfun(stderr.beta1.aghq[ITER,JTER,inds]),dp)
         mstderr.beta1.gva  <- round(locfun(stderr.beta1.gva[ITER,JTER,inds]),dp)
         
         rmse.beta1.pql  <- round(sqrt(locfun((beta1.pql[ITER,JTER,inds] - beta1.true)^2)),dp)
         rmse.beta1.lap  <- round(sqrt(locfun((beta1.lap[ITER,JTER,inds] - beta1.true)^2)),dp)
         rmse.beta1.aghq <- round(sqrt(locfun((beta1.aghq[ITER,JTER,inds] - beta1.true)^2)),dp)
         rmse.beta1.gva  <- round(sqrt(locfun((beta1.gva[ITER,JTER,inds] - beta1.true)^2)),dp)       
                  
         cat(mean.beta1.pql, "(",estderr.beta1.pql ,") & " ,"(",sd.beta1.pql ,") & ", mstderr.beta1.pql ," & ", rmse.beta1.pql, " \n")
         cat(mean.beta1.lap, "(",estderr.beta1.lap ,") & " ,"(",sd.beta1.lap ,") & ", mstderr.beta1.lap ," & ", rmse.beta1.lap, " \n")
         cat(mean.beta1.aghq,"(",estderr.beta1.aghq ,") & " ,"(",sd.beta1.aghq,") & ", mstderr.beta1.aghq," & ", rmse.beta1.aghq, " \n")
         cat(mean.beta1.gva, "(",estderr.beta1.gva ,") & " ,"(",sd.beta1.gva ,") & ", mstderr.beta1.gva ," & ", rmse.beta1.gva, " \n\n")  
         
         ######################################################################
            
         # Calculate and print summary statistic for fitted sigma values
                          
         cat("sigma \n")
         
         mean.sigma.pql  <- round(locfun(sigma.pql[ITER,JTER,inds]),dp)
         mean.sigma.lap  <- round(locfun(sigma.lap[ITER,JTER,inds]),dp)
         mean.sigma.aghq <- round(locfun(sigma.aghq[ITER,JTER,inds]),dp)
         mean.sigma.gva  <- round(locfun(sigma.gva[ITER,JTER,inds]),dp)
         
         estderr.sigma.pql  <- round(sprfun(sigma.pql[ITER,JTER,inds])/sqrt(length(inds)),dp)
         estderr.sigma.lap  <- round(sprfun(sigma.lap[ITER,JTER,inds])/sqrt(length(inds)),dp)
         estderr.sigma.aghq <- round(sprfun(sigma.aghq[ITER,JTER,inds]/sqrt(length(inds))),dp)
         estderr.sigma.gva  <- round(sprfun(sigma.gva[ITER,JTER,inds]/sqrt(length(inds))),dp)            
         
         sd.sigma.pql  <- round(sprfun(sigma.pql[ITER,JTER,inds]),dp)
         sd.sigma.lap  <- round(sprfun(sigma.lap[ITER,JTER,inds]),dp)
         sd.sigma.aghq <- round(sprfun(sigma.aghq[ITER,JTER,inds]),dp)
         sd.sigma.gva  <- round(sprfun(sigma.gva[ITER,JTER,inds]),dp)  
                 
         rmse.sigma.pql  <- round(sqrt(locfun((sigma.pql[ITER,JTER,inds] - sigma.true)^2)),dp)
         rmse.sigma.lap  <- round(sqrt(locfun((sigma.lap[ITER,JTER,inds] - sigma.true)^2)),dp)
         rmse.sigma.aghq <- round(sqrt(locfun((sigma.aghq[ITER,JTER,inds] - sigma.true)^2)),dp)
         rmse.sigma.gva  <- round(sqrt(locfun((sigma.gva[ITER,JTER,inds] - sigma.true)^2)),dp)  
         
         cat(mean.sigma.pql ,"(",estderr.sigma.pql  ,") &  (",sd.sigma.pql ,") &  ", rmse.sigma.pql, " \n")
         cat(mean.sigma.lap ,"(",estderr.sigma.lap  ,") &  (",sd.sigma.lap ,") &  ", rmse.sigma.lap, " \n")
         cat(mean.sigma.aghq,"(",estderr.sigma.aghq ,") &  (",sd.sigma.aghq,") &  ", rmse.sigma.aghq, " \n")
         cat(mean.sigma.gva ,"(",estderr.sigma.gva  ,") &  (",sd.sigma.gva ,") &  ", rmse.sigma.gva, " \n\n")  
         
         ######################################################################
         
         # Display mean times         
                  
         cat("times - ", round(mean(times.pql[ITER,JTER,inds]),dp), 
                         round(mean(times.lap[ITER,JTER,inds]),dp), 
                         round(mean(times.aghq[ITER,JTER,inds]),dp), 
                         round(mean(times.gva[ITER,JTER,inds]),dp),"\n")   
                         
         ######################################################################
         
         # Display the mean Root Mean Square Error (RMSE) for predicted random effects                         
                         
         dp <- 3 
                  
         inds2 <- is.included[is.included<TRIALS2]                  
                           
         cat("rmse.vmu - ",round(mean(rmse.vmu.lap[ITER,JTER,inds2]),dp), 
                           round(mean(rmse.vmu.gva[ITER,JTER,inds2]),dp),"\n")                            
          
         
         cat("\n\n")        
         
         ###################################################################### 
     }
}

###############################################################################

