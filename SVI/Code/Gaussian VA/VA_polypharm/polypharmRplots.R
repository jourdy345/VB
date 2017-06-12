library(rstan)


setwd("C:\\Users\\statsll\\Dropbox\\Gaussian VA revised\\VA_polypharm")


#####################
# polypharmacy data #
#####################

df <-  read.table("data_polypharm.txt", header=TRUE)
ID <- unique(df$ID)
n <- length(ID)                      # no. of subjects
vni <- rep(0,n)
X <- NULL
Z <- NULL
y <- NULL
N <- dim(df)[1]
AGEFP1 <- log(df$AGE/10)
gender_M <- df$GENDER
MHV4_1 <- rep(0,N)
MHV4_2 <- rep(0,N)
MHV4_3 <- rep(0,N)
RACE2_1 <- rep(0,N)
INPTMHV2_1 <- rep(0,N)
MHV4_1[df$MHV4==1] <- 1 
MHV4_2[df$MHV4==2] <- 1
MHV4_3[df$MHV4==3] <- 1
RACE2_1[df$RACE>0] <- 1
INPTMHV2_1[df$INPTMHV3>0] <- 1

for (i in 1:n){
    rows <- which(df$ID == ID[i])
    vni[i] <- length(rows)
    y[[i]] <- (df$POLYPHARMACY)[rows]
    Z[[i]] <- cbind(rep(1,vni[i]))
    X[[i]] <- cbind(rep(1,vni[i]), gender_M[rows], RACE2_1[rows], AGEFP1[rows],
              MHV4_1[rows], MHV4_2[rows],MHV4_3[rows],INPTMHV2_1[rows])
}
sum(vni)                                    # total no. of observations
n <- length(y)                              # no. of subjects
k <- dim(X[[1]])[2]                         # length of beta
p <- dim(Z[[1]])[2]                         # length of b_i 
d <- n*p + k + p*(p+1)/2
yall <- unlist(y)
Xall <- matrix(unlist(lapply(X,t)), ncol=k, byrow=TRUE)
Zall <- unlist(Z)
N <- sum(vni)     
startindex <- c(0, cumsum(vni)[1:(n-1)]) + 1
endindex <- cumsum(vni)
data <- list(n=n, N=N, k=k, y=yall, X=Xall, Z=Zall, startindex=startindex, endindex=endindex)
 


# stan #
fit <- stan(file = 'model_polypharm.stan', pars = c("beta", "zeta"),
            data = data, iter = 50000, chains = 1, thin=5)
print(fit, digits = 3)
la <- extract(fit, permuted = TRUE) # return a list of arrays 
mcmc <- cbind(la$beta, la$zeta)
par(mfrow=c(5,2))
par(mar=c(1,2,1,1))
for (i in 1:9){ plot(mcmc[,i],type='l') }
par(mfrow=c(2,5))
for (i in 1:9){ plot(density(mcmc[,i]))}
mcmc_mean <- apply(mcmc,2,mean)
mcmc_sd <- apply(mcmc,2,sd)
# write.table(mcmc, file = "mcmc_polypharm.txt", row.names = FALSE, col.names = FALSE)

#  Elapsed Time: 1305.61 seconds (Warm-up)
#                1350.82 seconds (Sampling)
#                2656.44 seconds (Total)
             mean se_mean     sd      2.5%       25%       50%       75%     97.5% n_eff Rhat
beta[1]    -4.258   0.007  0.393    -5.077    -4.515    -4.253    -3.989    -3.520  3038    1
beta[2]     0.757   0.006  0.342     0.092     0.532     0.758     0.982     1.431  3205    1
beta[3]    -0.693   0.007  0.380    -1.442    -0.945    -0.688    -0.445     0.041  3338    1
beta[4]     2.619   0.005  0.303     2.036     2.411     2.617     2.824     3.217  4518    1
beta[5]     0.330   0.005  0.287    -0.224     0.134     0.324     0.523     0.900  3972    1
beta[6]     1.209   0.005  0.298     0.634     1.007     1.203     1.408     1.819  3792    1
beta[7]     1.749   0.005  0.296     1.165     1.546     1.750     1.944     2.345  4107    1
beta[8]     0.894   0.004  0.259     0.388     0.718     0.892     1.072     1.394  4994    1
zeta        0.916   0.002  0.068     0.785     0.870     0.916     0.961     1.052  1883    1
lp__    -1736.845   0.577 25.200 -1788.016 -1753.339 -1736.398 -1719.928 -1688.845  1906    1





# stan (compare mcmc, mean-field, fullrank)
model_polypharm <- stan_model(file='model_polypharm.stan')

# mean-field #
system.time(fit.mf <- vb(model_polypharm, pars=c("beta", "zeta"),
            data=data, algorithm="meanfield", output_samples=2000,seed=123))
print(fit.mf)
ex.mf <- extract(fit.mf, permuted = TRUE)
vb.mf <- cbind(ex.mf$beta, ex.mf$zeta)
# save(vb.mf, file="stan_polypharm_vb.mf")
# Inference for Stan model: model_polypharm.
1 chains, each with iter=2000; warmup=0; thin=1; 
post-warmup draws per chain=2000, total post-warmup draws=2000.

   user  system elapsed 
   6.79    0.32    7.33 

         mean   sd  2.5%   25%   50%   75% 97.5%
beta[1] -4.01 0.05 -4.11 -4.04 -4.01 -3.98 -3.91
beta[2]  0.37 0.06  0.25  0.33  0.37  0.42  0.50
beta[3] -0.85 0.14 -1.13 -0.95 -0.86 -0.76 -0.57
beta[4]  2.66 0.20  2.27  2.52  2.66  2.80  3.04
beta[5]  0.28 0.14  0.02  0.19  0.29  0.38  0.56
beta[6]  1.14 0.10  0.95  1.07  1.14  1.20  1.33
beta[7]  1.70 0.10  1.51  1.64  1.70  1.77  1.90
beta[8]  0.95 0.23  0.51  0.79  0.94  1.10  1.41
zeta     0.92 0.03  0.86  0.90  0.92  0.94  0.98
lp__     0.00 0.00  0.00  0.00  0.00  0.00  0.00

Approximate samples were drawn using VB(meanfield) at Fri Nov 11 09:47:45 2016.
We recommend genuine 'sampling' from the posterior distribution for final inferences!



# full-rank
system.time(fit.fr <- vb(model_polypharm, pars = c("beta", "zeta"),
            data = data, algorithm="fullrank", output_samples=2000,seed=123))
print(fit.fr)
plot(fit.fr)
ex.fr <- extract(fit.fr, permuted = TRUE)
vb.fr <- cbind(ex.fr$beta, ex.fr$zeta)

# save(vb.fr, file="stan_polypharm_vb.fr")
# Inference for Stan model: model_polypharm.
1 chains, each with iter=2000; warmup=0; thin=1; 
post-warmup draws per chain=2000, total post-warmup draws=2000.

   user  system elapsed 
  48.52   25.35   74.29 

         mean   sd  2.5%   25%   50%   75% 97.5%
beta[1] -2.21 0.45 -3.14 -2.50 -2.21 -1.92 -1.30
beta[2] -0.27 0.37 -1.01 -0.50 -0.27 -0.01  0.44
beta[3] -1.10 0.31 -1.70 -1.30 -1.09 -0.90 -0.49
beta[4]  2.21 0.28  1.65  2.03  2.21  2.41  2.74
beta[5] -0.68 0.35 -1.36 -0.91 -0.69 -0.45  0.01
beta[6]  0.23 0.33 -0.41  0.01  0.23  0.45  0.89
beta[7]  0.87 0.33  0.22  0.65  0.86  1.10  1.54
beta[8]  1.07 0.28  0.54  0.89  1.06  1.26  1.63
zeta     0.91 0.05  0.81  0.88  0.91  0.95  1.01
lp__     0.00 0.00  0.00  0.00  0.00  0.00  0.00

Approximate samples were drawn using VB(fullrank) at Fri Nov 11 09:52:23 2016.
We recommend genuine 'sampling' from the posterior distribution for final inferences!




load("stan_polypharm_vb.mf")
load("stan_polypharm_vb.fr")
mcmc <- read.table("mcmc_polypharm.txt", header=F)
VA2ameansd <- read.table("VA2ameansd_intercept.txt", header=F)
VA1diagbmeansd <- read.table("VA1diagbmeansd_intercept.txt", header=F)
VA1bmeansd <- read.table("VA1bmeansd_intercept.txt", header=F)


labels <- c(expression(beta[0]), expression(beta[Gender]), expression(beta[Race]),
expression(beta[Age]), expression(beta[MHV4[1]]), expression(beta[MHV4[2]]), 
expression(beta[MHV4[3]]),expression(beta[INPTMHV2[1]]),expression(zeta[1]))
ymax <- c(2, 2, 2.9, 2, 2.9, 4, 4.5, 1.9, 13)

pdf("C:\\Users\\statsll\\Dropbox\\Gaussian VA revised\\S&C_Manuscript\\images\\polypharm_intercept.pdf", height=4, width=10)
par(mar=c(3.5,2.3,1,1))
par(mfrow=c(2,5))
par(mgp=c(2.7,0.7,0))
for (i in 1:9){
 temp <- density(mcmc[,i])
 plot(temp, main="", lwd=1.5, 
       ylab="", xlab=labels[i], ylim=c(0,ymax[i]), cex.lab=1.8, cex.axis=1.1)
 points(density(vb.fr[,i]),type="l",col="green", lwd=1.5, lty=2)
 points(density(vb.mf[,i]),type="l",col="blue", lwd=1.5, lty=2)
 lines(temp$x, dnorm(temp$x,mean=VA1diagbmeansd[i,1],sd=VA1diagbmeansd[i,2]), 
      type='l', lty=1,col="blue", lwd=1.5)
 lines(temp$x, dnorm(temp$x,mean=VA1bmeansd[i,1],sd=VA1bmeansd[i,2]), 
      type='l', lty=1,col="green", lwd=1.5)
 lines(temp$x, dnorm(temp$x,mean=VA2ameansd[i,1],sd=VA2ameansd[i,2]),
      type='l', lty=1,col="red", lwd=1.5)
}
dev.off()





