library(rstan)

setwd("C:\\Users\\statsll\\Dropbox\\Gaussian VA revised\\VA_toenail")

################
# toenail data #
################

df <-  read.table("data_toenail.txt", header=TRUE)
ID <- unique(df$ID)
n <- length(ID)                      # no. of subjects
vni <- rep(0,n)
X <- NULL
Z <- NULL
y <- NULL
for (i in 1:n){
    rows <- which(df$ID == ID[i])
    vni[i] <- length(rows)
    y[[i]] <- (df$Response)[rows]
    Z[[i]] <- cbind(rep(1,vni[i]))
    X[[i]] <- cbind(rep(1,vni[i]), df$Treatment[rows], df$Month[rows], df$Treatment[rows]*df$Month[rows])
}
sum(vni)                                    # total no. of observations
n <- length(y)                              # no. of subjects
k <- dim(X[[1]])[2]                         # length of beta
p <- dim(Z[[1]])[2]                         # length of b_i 
d <- n*p + k + p*(p+1)/2
yall <- unlist(y)
Xall <- matrix(unlist(lapply(X,t)), ncol=k, byrow=TRUE)
Zall <- unlist(Z)
labels <- c('intercept','trt','time','trt x time','zeta')  
N <- sum(vni)
startindex <- c(0, cumsum(vni)[1:(n-1)]) + 1
endindex <- cumsum(vni)
data <- list(n=n, N=N, k=k, y=yall, X=Xall, Z=Zall, startindex=startindex, endindex=endindex)



########
# stan #
########

fit <- stan(file = 'model_toenail.stan', pars = c("beta", "zeta"),
            data = data, iter = 50000, chains = 1, thin=5)
print(fit, digits = 3)
la <- extract(fit, permuted = TRUE) # return a list of arrays 
mcmc <- cbind(la$beta, la$zeta)

par(mfrow=c(5,1))
par(mar=c(1,2,1,1))
for (i in 1:5){ plot(mcmc[,i],type='l') }

par(mfrow=c(1,5))
for (i in 1:5){ plot(density(mcmc[,i]))}
mcmc_mean <- apply(mcmc,2,mean)
mcmc_sd <- apply(mcmc,2,sd)
#write.table(mcmc, file="mcmc_toenail.txt", row.names=FALSE, col.names=FALSE)


#  Elapsed Time: 441.627 seconds (Warm-up)
#                664.752 seconds (Sampling)
#                1106.38 seconds (Total)
            mean se_mean     sd     2.5%      25%      50%      75%    97.5% n_eff  Rhat
beta[1]   -1.662   0.008  0.445   -2.577   -1.956   -1.649   -1.350   -0.831  3054 1.001
beta[2]   -0.160   0.009  0.604   -1.353   -0.562   -0.156    0.255    1.014  4284 1.000
beta[3]   -0.398   0.001  0.045   -0.496   -0.427   -0.396   -0.367   -0.314  4169 1.000
beta[4]   -0.140   0.001  0.069   -0.274   -0.186   -0.138   -0.094   -0.007  4792 1.000
zeta       1.414   0.002  0.094    1.229    1.350    1.413    1.477    1.601  2579 1.000
lp__    -952.875   0.455 21.512 -995.420 -967.118 -952.235 -938.114 -912.424  2233 1.000



# stan (compare mcmc, mean-field, fullrank)
model_toenail <- stan_model(file='model_toenail.stan')

# mean-field #
system.time(fit.mf <- vb(model_toenail, pars=c("beta", "zeta"),
            data=data, algorithm="meanfield", output_samples=2000,seed=123))
print(fit.mf)
ex.mf <- extract(fit.mf, permuted = TRUE)
vb.mf <- cbind(ex.mf$beta, ex.mf$zeta)
# save(vb.mf, file="stan_toenail_vb.mf")
# Inference for Stan model: model_toenail.
# 1 chains, each with iter=2000; warmup=0; thin=1; 
# post-warmup draws per chain=2000, total post-warmup draws=2000.

   user  system elapsed 
   3.98    0.12    4.14 


         mean   sd  2.5%   25%   50%   75% 97.5%
beta[1] -1.32 0.09 -1.50 -1.39 -1.32 -1.26 -1.14
beta[2] -0.30 0.13 -0.55 -0.38 -0.30 -0.21 -0.03
beta[3] -0.35 0.03 -0.40 -0.37 -0.35 -0.33 -0.30
beta[4] -0.14 0.03 -0.20 -0.16 -0.14 -0.12 -0.08
zeta     1.29 0.04  1.21  1.26  1.29  1.32  1.38
lp__     0.00 0.00  0.00  0.00  0.00  0.00  0.00

Approximate samples were drawn using VB(meanfield) at Fri Nov 11 09:27:40 2016.


# full-rank
system.time(fit.fr <- vb(model_toenail, pars = c("beta", "zeta"),
            data = data, algorithm="fullrank", output_samples=2000,seed=123))
print(fit.fr)
plot(fit.fr)
ex.fr <- extract(fit.fr, permuted = TRUE)
vb.fr <- cbind(ex.fr$beta, ex.fr$zeta)

# save(vb.fr, file="stan_toenail_vb.fr")
# Inference for Stan model: model_toenail.
# 1 chains, each with iter=2000; warmup=0; thin=1; 
# post-warmup draws per chain=2000, total post-warmup draws=2000.

   user  system elapsed 
  12.03    3.20   15.26 

         mean   sd  2.5%   25%   50%   75% 97.5%
beta[1] -0.39 0.36 -1.08 -0.64 -0.39 -0.14  0.31
beta[2] -0.09 0.40 -0.85 -0.36 -0.10  0.16  0.70
beta[3] -0.32 0.05 -0.42 -0.35 -0.32 -0.29 -0.23
beta[4] -0.08 0.07 -0.21 -0.13 -0.08 -0.04  0.04
zeta     0.93 0.08  0.78  0.88  0.93  0.98  1.08
lp__     0.00 0.00  0.00  0.00  0.00  0.00  0.00

Approximate samples were drawn using VB(fullrank) at Fri Nov 11 09:28:58 2016.




load("stan_toenail_vb.mf")
load("stan_toenail_vb.fr")
mcmc <- read.table("mcmc_toenail.txt", header=F)
VA2ameansd <- read.table("VA2ameansd_intercept.txt", header=F)
VA1diagbmeansd <- read.table("VA1diagbmeansd_intercept.txt", header=F)
VA1bmeansd <- read.table("VA1bmeansd_intercept.txt", header=F)


labels <- c(expression(beta[0]), expression(beta[Trt]), expression(beta[t]),
expression(beta[Trtxt]),expression(zeta[1]))


# pdf("C:\\Users\\statsll\\Dropbox\\Gaussian VA revised\\S&C_Manuscript\\images\\toenail_intercept.pdf", height=1.8, width=8.5)
ymax <- c(4.2, 2.7, 18, 13, 9.7)
par(mar=c(3.5,2.1,1,0.9))
par(mfrow=c(1,5))
par(mgp=c(2.7,0.7,0))
for (i in 1:5){
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





ymin = c(-0.2, -0.1, -1, -0.2, -0.1, -0.35, -0.7, -0.1, -0.3)
ymax = c(0.4, 0.9, 0.2, 0.4, 0.6, 0, 0.5, 0.07, 0.5)
pdf("C:\\Users\\Linda\\Desktop\\Gaussian VA\\script\\toenail_slope_convergence.pdf", height=3.5, width=10)
par(mfcol=c(2,5))
par(mar=c(3.5,2.3,2,1))
par(mgp=c(2.2,0.7,0))
for (i in 1:5){
 plot(VA2apar[,i], type="l", lwd=1.5, ,col="green", main = labels[i], #ylim=c(ymin[i], ymax[i]),
       ylab="", xlab= "iterations (x 1000)", cex.main=1.5, cex.axis=1.1, cex.lab=1.1)
 plot(VA1bpar[,i], type='l', lty=1,col="red", lwd=1.5)
}
dev.off(




gclosed2a <- read.table("VA2agclosed.txt", header=F)
gapprox2a <- read.table("VA2agapprox.txt", header=F)


pdf("C:\\Users\\statsll\\Dropbox\\Gaussian VA revised\\S&C_Manuscript\\images\\toenail_Alg2gradtest.pdf", height=2, width=8.5)
par(mar=c(3.3,2.3,0.8,0.8))
par(mfrow=c(1,5))
par(mgp=c(2.5,0.7,0))
for (i in 1:5){
 plot(1:1000, gclosed2a[,i], type='l', lwd=1, ylab="", xlab=labels[i], cex.lab=1.7, cex.axis=1)
 lines(1:1000, gapprox2a[,i], type='l', lty=1,col="red", lwd=1)
}
dev.off()




