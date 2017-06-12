library(rstan)

setwd("C:\\Users\\statsll\\Dropbox\\Gaussian VA revised\\VA_epilepsy")



####################
# Epilepsy dataset #
####################

df <-  read.table("data_epilepsy.txt", header=TRUE)
n <- dim(df)[1]                      # no. of subjects
vni <- rep(4,n)
X <- NULL
Z <- NULL
y <- NULL
lb4 <- log(df$Base4)
trt <- c(rep(0,28),rep(1,31))
lage <- log(df$Age)
clage <- scale(lage,scale=FALSE)
V4 <- c(0,0,0,1)
visit <- c(-3,-1,1,3)/10
lb4trt <- lb4*trt
N <- sum(vni)
startindex <- c(0, cumsum(vni)[1:(n-1)]) + 1
endindex <- cumsum(vni)




################
# random slope #
################
for (i in 1:n){
    y[[i]] <- c(df$Y1[i], df$Y2[i], df$Y3[i], df$Y4[i])
    Z[[i]] <- cbind(rep(1,vni[i]), visit)
    X[[i]] <- cbind(rep(1,vni[i]), rep(lb4[i],4), rep(trt[i],4), rep(lb4trt[i],4), rep(clage[i],4), visit)
}
sum(vni)                                    # total no. of observations
n <- length(y)                              # no. of subjects
k <- dim(X[[1]])[2]                         # length of beta
p <- dim(Z[[1]])[2]                         # length of b_i 
d <- n*p + k + p*(p+1)/2
yall <- unlist(y)
Xall <- matrix(unlist(lapply(X,t)), ncol=k, byrow=TRUE)
Zall <- matrix(unlist(lapply(Z,t)), ncol=p, byrow=TRUE)
pzeros <- rep(0,p)
zetalength <- p*(p+1)/2
data <- list(n=n, N=N, k=k, p=p, y=yall, X=Xall, Z=Zall, pzeros=pzeros,
             startindex=startindex, endindex=endindex, zetalength=zetalength)


# stan (mcmc) #
fit <- stan(file = 'model_epilepsy_slope.stan', pars = c("beta", "zeta"),
            data = data, iter = 50000, chains = 1,thin=5)


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
# write.table(mcmc, file = "mcmc_epilepsy_slope.txt", row.names = FALSE, col.names = FALSE)

#  Elapsed Time: 768.91 seconds (Warm-up)
#                1310.11 seconds (Sampling)
#                2079.02 seconds (Total)
            mean se_mean     sd     2.5%      25%      50%      75%    97.5% n_eff  Rhat
beta[1]    0.205   0.004  0.269   -0.321    0.025    0.201    0.382    0.738  3953 1.000
beta[2]    0.886   0.002  0.137    0.618    0.795    0.888    0.978    1.152  3827 1.000
beta[3]   -0.942   0.006  0.422   -1.762   -1.221   -0.941   -0.659   -0.116  4288 1.000
beta[4]    0.343   0.003  0.215   -0.084    0.206    0.343    0.488    0.760  4018 1.001
beta[5]    0.476   0.006  0.377   -0.271    0.223    0.475    0.727    1.214  3972 1.000
beta[6]   -0.269   0.002  0.169   -0.608   -0.382   -0.272   -0.158    0.065  4642 1.000
zeta[1]   -0.612   0.002  0.121   -0.842   -0.697   -0.614   -0.531   -0.369  4768 1.000
zeta[2]    0.004   0.003  0.185   -0.359   -0.118    0.002    0.129    0.367  4085 1.000
zeta[3]   -0.322   0.010  0.273   -0.825   -0.459   -0.302   -0.157    0.103   705 1.000
lp__    3222.456   0.508 13.994 3198.849 3213.952 3221.976 3229.863 3248.281   759 1.000



# stan (compare mcmc, mean-field, fullrank)
model_epilepsy_slope <- stan_model(file='model_epilepsy_slope.stan')

# mean-field #
system.time(fit.mf <- vb(model_epilepsy_slope, 
                      pars=c("beta", "zeta"),
                      data=data, algorithm="meanfield", 
                      output_samples=2000,seed=123) )
   user  system elapsed 
   4.15    0.22    4.38 

print(fit.mf)
ex.mf <- extract(fit.mf, permuted = TRUE)
vb.mf <- cbind(ex.mf$beta, ex.mf$zeta)
# save(vb.mf, file="stan_slope_vb.mf")
# Inference for Stan model: model_epilepsy_slope.
# 1 chains, each with iter=2000; warmup=0; thin=1; 
# post-warmup draws per chain=2000, total post-warmup draws=2000.

         mean   sd  2.5%   25%   50%   75% 97.5%
beta[1]  0.89 0.06  0.76  0.85  0.89  0.94  1.02
beta[2]  0.58 0.03  0.52  0.56  0.58  0.60  0.64
beta[3] -0.84 0.08 -1.00 -0.89 -0.83 -0.78 -0.67
beta[4]  0.73 0.04  0.65  0.70  0.73  0.75  0.80
beta[5]  1.18 0.15  0.90  1.08  1.19  1.29  1.48
beta[6] -0.06 0.12 -0.29 -0.14 -0.06  0.02  0.16
zeta[1] -0.04 0.10 -0.23 -0.10 -0.04  0.03  0.15
zeta[2]  0.24 0.13 -0.01  0.15  0.24  0.31  0.49
zeta[3] -0.11 0.10 -0.30 -0.18 -0.11 -0.05  0.08
lp__     0.00 0.00  0.00  0.00  0.00  0.00  0.00

Approximate samples were drawn using VB(meanfield) at Fri Nov 11 08:36:49 2016.
We recommend genuine 'sampling' from the posterior distribution for final inferences!


# full-rank
system.time(fit.fr <- vb(model_epilepsy_slope, 
                      pars = c("beta", "zeta"),
                      data = data, algorithm="fullrank", 
                      output_samples=2000,seed=123) )
   user  system elapsed 
   3.77    0.16    3.92 

print(fit.fr)
plot(fit.fr)
ex.fr <- extract(fit.fr, permuted = TRUE)
vb.fr <- cbind(ex.fr$beta, ex.fr$zeta)

# save(vb.fr, file="stan_slope_vb.fr")
# Inference for Stan model: model_epilepsy_slope.
# 1 chains, each with iter=2000; warmup=0; thin=1; 
# post-warmup draws per chain=2000, total post-warmup draws=2000.

         mean   sd  2.5%   25%   50%   75% 97.5%
beta[1]  0.80 0.73 -0.67  0.31  0.81  1.32  2.21
beta[2]  0.59 0.46 -0.26  0.28  0.60  0.90  1.48
beta[3] -0.95 0.78 -2.59 -1.47 -0.95 -0.43  0.59
beta[4]  0.73 0.42 -0.11  0.45  0.74  1.02  1.52
beta[5]  0.79 0.81 -0.72  0.24  0.76  1.34  2.40
beta[6]  0.03 0.36 -0.65 -0.22  0.04  0.26  0.74
zeta[1]  0.19 0.21 -0.23  0.06  0.18  0.33  0.59
zeta[2]  0.17 0.29 -0.42 -0.03  0.18  0.36  0.74
zeta[3]  0.19 0.13 -0.07  0.10  0.19  0.28  0.44
lp__     0.00 0.00  0.00  0.00  0.00  0.00  0.00

Approximate samples were drawn using VB(fullrank) at Fri Nov 11 08:38:04 2016.
We recommend genuine 'sampling' from the posterior distribution for final inferences!



# Plots
load("stan_slope_vb.mf")
load("stan_slope_vb.fr")
mcmc <- read.table(file="mcmc_epilepsy_slope.txt")
VA2ameansd <- read.table("VA2ameansd_slope.txt", header=F)
VA1diagbmeansd <- read.table("VA1diagbmeansd_slope.txt", header=F)
VA1bmeansd <- read.table("VA1bmeansd_slope.txt", header=F)

labels <- c(expression(beta[0]), expression(beta[Base]), expression(beta[Trt]),
expression(beta[BasexTrt]), expression(beta[Age]), expression(beta[Visit]), 
expression(zeta[1]),expression(zeta[2]),expression(zeta[3]),"Lower bound")


# pdf("C:\\Users\\statsll\\Dropbox\\Gaussian VA revised\\S&C_Manuscript\\images\\epilepsy_slope.pdf",
 height=6, width=9)
par(mar=c(3.5,2.3,1,1))
par(mfrow=c(3,3))
par(mgp=c(2.7,0.7,0))
ymax <- c(2, 3.5, 1.5, 2.5, 1.5, 3.9, 4.2, 3.6, 4.2)
xmin=c(-1.1,0.3,-2.5,-0.5,-1,-0.8,-1.1,-0.7,-1.1)
xmax=c(1.3,1.4,0.7,1.2,2,0.4,-0.1,0.7,0.6)
for (i in 1:9){
 temp <- density(mcmc[,i])
 plot(temp, main="", ylim=c(0,ymax[i]), xlim=c(xmin[i],xmax[i]),
       ylab="", xlab=labels[i], cex.lab=1.7, cex.axis=1.1, lwd=1.5)
 points(density(vb.fr[,i]),type="l",col="green", lwd=1.5, lty=2)
 points(density(vb.mf[,i]),type="l",col="blue", lwd=1.5, lty=2)
 lines(temp$x, dnorm(temp$x,mean=VA1diagbmeansd[i,1],sd=VA1diagbmeansd[i,2]), 
      type='l',col="blue", lwd=1.5)
 lines(temp$x, dnorm(temp$x,mean=VA1bmeansd[i,1],sd=VA1bmeansd[i,2]), 
      type='l',col="green", lwd=1.5)
 lines(temp$x, dnorm(temp$x,mean=VA2ameansd[i,1],sd=VA2ameansd[i,2]),
      type='l',col="red", lwd=1.5)
}
dev.off()


VA2apar <- read.table("VA2apar_slope.txt", header=F)
VA1bpar <- read.table("VA1bpar_slope.txt", header=F)

# pdf("C:\\Users\\statsll\\Dropbox\\Gaussian VA revised\\S&C_Manuscript\\images\\epilepsy_slope_convergence.pdf",
height=4, width=9.5)
ymin = c(-0.2, -0.1, -1, -0.2, -0.1, -0.35, -0.7, -0.1, -0.3,0)
ymax = c(0.4, 0.9, 0.2, 0.4, 0.6, 0, 0.5, 0.07, 0.5,3250)
par(mfrow=c(2,5))
par(mar=c(3.5,2.3,2,1))
par(mgp=c(2.2,0.7,0))
j <- 0
for (i in c(1:9,19)){
 j <- j+1
 plot(VA2apar[,i], type="o", lwd=1, ,col="green", main = labels[j], ylim=c(ymin[j], ymax[j]),
       ylab="", xlab= "iterations (x 2500)", cex.main=1.5, cex.axis=1.1, cex.lab=1.1,pch=20)
 lines(VA1bpar[,i], type='o', lty=1,col="red", lwd=1,pch=20)
}
dev.off()











####################
# random intercept #
####################
for (i in 1:n){
    y[[i]] <- c(df$Y1[i], df$Y2[i], df$Y3[i], df$Y4[i])
    Z[[i]] <- cbind(rep(1,vni[i]))
    X[[i]] <- cbind(rep(1,vni[i]), rep(lb4[i],4), rep(trt[i],4), rep(lb4trt[i],4), rep(clage[i],4), V4)
}
sum(vni)                                    # total no. of observations
n <- length(y)                              # no. of subjects
k <- dim(X[[1]])[2]                         # length of beta
p <- dim(Z[[1]])[2]                         # length of b_i 
d <- n*p + k + p*(p+1)/2
yall <- unlist(y)
Xall <- matrix(unlist(lapply(X,t)), ncol=k, byrow=TRUE)
Zall <- unlist(Z)
data <- list(n=n, N=N, k=k, y=yall, X=Xall, Z=Zall,
             startindex=startindex, endindex=endindex, zetalength=zetalength)


# stan (mcmc)#
fit <- stan(file = 'model_epilepsy_intercept.stan', pars = c("beta", "zeta"),
            data = data, iter = 50000, chains = 1, thin=5)

print(fit, digits = 3)
la <- extract(fit, permuted = TRUE) # return a list of arrays 
mcmc <- cbind(la$beta, la$zeta)
par(mfrow=c(4,2))
par(mar=c(1,2,1,1))
for (i in 1:7){ plot(mcmc[,i],type='l') }
par(mfrow=c(2,4))
for (i in 1:7){ plot(density(mcmc[,i]))}
mcmc_mean <- apply(mcmc,2,mean)
mcmc_sd <- apply(mcmc,2,sd)
# write.table(mcmc, file = "mcmc_epilepsy_intercept.txt", row.names = FALSE, col.names = FALSE)

#  Elapsed Time: 184.433 seconds (Warm-up)
#                341.807 seconds (Sampling)
#                526.24 seconds (Total)
            mean se_mean    sd     2.5%      25%      50%      75%    97.5% n_eff Rhat
beta[1]    0.269   0.004 0.273   -0.263    0.084    0.270    0.452    0.796  4483    1
beta[2]    0.883   0.002 0.140    0.609    0.791    0.883    0.974    1.161  4315    1
beta[3]   -0.940   0.006 0.426   -1.776   -1.222   -0.936   -0.655   -0.113  4571    1
beta[4]    0.339   0.003 0.219   -0.092    0.194    0.339    0.484    0.774  4384    1
beta[5]    0.479   0.005 0.371   -0.252    0.236    0.482    0.728    1.204  4638    1
beta[6]   -0.161   0.001 0.054   -0.268   -0.197   -0.160   -0.124   -0.055  4916    1
zeta      -0.621   0.002 0.120   -0.849   -0.704   -0.623   -0.538   -0.380  4716    1
lp__    3208.535   0.092 6.168 3195.651 3204.664 3208.827 3212.886 3219.654  4500    1




# stan (compare mcmc, mean-field, fullrank)
model_epilepsy_intercept <- stan_model(file='model_epilepsy_intercept.stan')

# mean-field #
system.time(fit.mf <- vb(model_epilepsy_intercept, pars=c("beta", "zeta"),
            data=data, algorithm="meanfield", output_samples=2000,seed=123))
  user  system elapsed 
   1.37    0.14    1.54 
 

print(fit.mf)
ex.mf <- extract(fit.mf, permuted = TRUE)
vb.mf <- cbind(ex.mf$beta, ex.mf$zeta)
# save(vb.mf, file="stan_intercept_vb.mf")
# Inference for Stan model: model_epilepsy_intercept.
# 1 chains, each with iter=2000; warmup=0; thin=1; 
# post-warmup draws per chain=2000, total post-warmup draws=2000.

         mean   sd  2.5%   25%   50%   75% 97.5%
beta[1]  0.58 0.02  0.54  0.57  0.58  0.59  0.62
beta[2]  0.72 0.01  0.70  0.71  0.72  0.72  0.74
beta[3] -1.43 0.03 -1.49 -1.45 -1.43 -1.41 -1.37
beta[4]  0.61 0.01  0.58  0.60  0.61  0.62  0.63
beta[5]  0.49 0.12  0.26  0.41  0.49  0.58  0.72
beta[6] -0.13 0.06 -0.24 -0.17 -0.13 -0.09 -0.02
zeta    -0.67 0.08 -0.84 -0.73 -0.67 -0.62 -0.50
lp__     0.00 0.00  0.00  0.00  0.00  0.00  0.00

Approximate samples were drawn using VB(meanfield) at Fri Nov 11 08:57:11 2016.
We recommend genuine 'sampling' from the posterior distribution for final inferences!



# full-rank
system.time(fit.fr <- vb(model_epilepsy_intercept, pars = c("beta", "zeta"),
            data = data, algorithm="fullrank", output_samples=2000,seed=123))
   user  system elapsed 
   1.21    0.19    1.37 

ex.fr <- extract(fit.fr, permuted = TRUE)
vb.fr <- cbind(ex.fr$beta, ex.fr$zeta)
print(fit.fr)
plot(fit.fr)
# save(vb.fr, file="stan_intercept_vb.fr")
# Inference for Stan model: model_epilepsy_intercept.
# 1 chains, each with iter=2000; warmup=0; thin=1; 
# post-warmup draws per chain=2000, total post-warmup draws=2000.

         mean   sd  2.5%   25%   50%   75% 97.5%
beta[1]  1.03 0.79 -0.52  0.48  1.03  1.57  2.54
beta[2]  0.65 0.44 -0.24  0.35  0.65  0.94  1.51
beta[3] -0.80 0.81 -2.49 -1.33 -0.80 -0.26  0.73
beta[4]  0.65 0.41 -0.19  0.38  0.65  0.92  1.41
beta[5]  1.25 1.08 -0.90  0.51  1.22  1.99  3.42
beta[6] -0.16 0.06 -0.27 -0.20 -0.16 -0.12 -0.05
zeta     0.28 0.27 -0.25  0.09  0.27  0.46  0.81
lp__     0.00 0.00  0.00  0.00  0.00  0.00  0.00

Approximate samples were drawn using VB(fullrank) at Fri Nov 11 08:57:51 2016.
We recommend genuine 'sampling' from the posterior distribution for final inferences!



# epilepsy_plots #

####################
# random intercept #
####################
load("stan_intercept_vb.mf")
load("stan_intercept_vb.fr")
VA2ameansd <- read.table("VA2ameansd_intercept.txt", header=F)
VA1diagbmeansd <- read.table("VA1diagbmeansd_intercept.txt", header=F)
VA1bmeansd <- read.table("VA1bmeansd_intercept.txt", header=F)
mcmc <- read.table("mcmc_epilepsy_intercept.txt", header=F)

labels <- c(expression(beta[0]), expression(beta[Base]), expression(beta[Trt]),
expression(beta[BasexTrt]), expression(beta[Age]), expression(beta[V4]), 
expression(zeta[1]))
ymax <- c(1.8, 3.7, 1.2, 2.4, 1.3, 8.5, 4.7)

# pdf("C:\\Users\\Li\\Dropbox\\Gaussian VA revised\\S&C_Manuscript\\images\\epilepsy_intercept.pdf", height=3.5, width=9)
par(mar=c(3.5,2.3,1,1))
par(mfrow=c(2,4))
par(mgp=c(2.7,0.7,0))
for (i in 1:7){
 temp <- density(mcmc[,i])
 plot(temp, main="", lwd=1.5, 
       ylab="", xlab=labels[i], ylim=c(0,ymax[i]), cex.lab=1.7, cex.axis=1.1)
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




