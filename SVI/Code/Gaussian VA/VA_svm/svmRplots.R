library(rstan)


setwd("C:\\Users\\statsll\\Dropbox\\Gaussian VA revised\\VA_svm")


#############################
# pounds exchange rate data #
#############################

rm(list=ls())


require(Ecdat)
data(Garch)
ex <- Garch$bp[which(Garch$date==811001):which(Garch$date==850628)]
T <- length(ex)
tex <- log(ex[2:T]) - log(ex[1:(T-1)]) 
y <- 100*(tex -mean(tex))
plot(y,type="l")
n <- length(y)
# write.table(y, file="exchange.txt",row.names=FALSE,col.names=FALSE)
data <- list(n=n, y=y)

# stan (mcmc) #
fit <- stan(file = 'model_svm2.stan', pars = c("alpha","lambda","psi","x"),
            data = data, iter = 50000, chains = 1, thin=5)
print(fit, digits = 3)
la <- extract(fit, permuted = TRUE) # return a list of arrays 
mcmc <- cbind(la$alpha, la$lambda, la$psi, la$x)
(mcmc_mean <- apply(mcmc,2,mean))
(mcmc_sd <- apply(mcmc,2,sd))
# save(mcmc, file="mcmc_exchange2")

par(mfrow=c(3,1))
par(mar=c(3,2,1,1))
for (i in 1:3){ plot(mcmc[,i],type='l') }

par(mfrow=c(1,3))
for (i in 1:3){ plot(density(mcmc[,i]))}

head(mcmc_mean)
head(mcmc_sd)

#  Elapsed Time: 1125.05 seconds (Warm-up)
#                1218.98 seconds (Sampling)
#                2344.02 seconds (Total)
[1] -1.8691120 -0.7091177  3.8685583  4.2824663  4.4688359  4.4712635
> head(mcmc_sd)
[1] 0.3194588 0.3792035 0.8918990 4.1454900 4.0580514 4.0100265




# stan (compare mcmc, mean-field, fullrank)
model_svm <- stan_model(file='model_svm2.stan')

# mean-field #
system.time(fit.mf <- vb(model_svm, pars = c("alpha","lambda","psi","x"),
            data=data, algorithm="meanfield", output_samples=2000,seed=123))


   user  system elapsed 
   9.40    0.34    9.76 

print(fit.mf)
ex.mf <- extract(fit.mf, permuted = TRUE)
vb.mf <-  cbind(ex.mf$alpha, ex.mf$lambda, ex.mf$psi, ex.mf$x)
# save(vb.mf, file="stan_pd_vb.mf")
# Inference for Stan model: model_svm2.
# 1 chains, each with iter=2000; warmup=0; thin=1; 
# post-warmup draws per chain=2000, total post-warmup draws=2000.

        mean   sd  2.5%   25%   50%   75% 97.5%
alpha  -0.42 0.06 -0.53 -0.46 -0.42 -0.38 -0.30
lambda -0.76 0.05 -0.86 -0.79 -0.76 -0.72 -0.66
psi    -3.89 1.06 -5.98 -4.61 -3.88 -3.17 -1.80

Approximate samples were drawn using VB(meanfield) at Fri Nov 11 09:47:45 2016.
We recommend genuine 'sampling' from the posterior distribution for final inferences!


# full-rank
system.time(fit.fr <- vb(model_svm, pars = c("alpha","lambda","psi","x"),
            data = data, algorithm="fullrank", output_samples=2000,seed=123))



print(fit.fr)
plot(fit.fr)
ex.fr <- extract(fit.fr, permuted = TRUE)
vb.fr <-  cbind(ex.mf$alpha, ex.mf$lambda, ex.mf$psi, ex.mf$x)

# This is Automatic Differentiation Variational Inference.

(EXPERIMENTAL ALGORITHM: expect frequent updates to the procedure.)

Gradient evaluation took 0 seconds
1000 iterations under these settings should take 0 seconds.
Adjust your expectations accordingly!

Begin eta adaptation.
Iteration:   1 / 250 [  0%]  (Adaptation)
Iteration:  50 / 250 [ 20%]  (Adaptation)
Iteration: 100 / 250 [ 40%]  (Adaptation)
Iteration: 150 / 250 [ 60%]  (Adaptation)
Iteration: 200 / 250 [ 80%]  (Adaptation)
Iteration: 250 / 250 [100%]  (Adaptation)
Error: stan::variational::advi::adapt_eta: All proposed step-sizes failed. Your model may be either severely ill-conditioned or misspecified.
Timing stopped at: 8.79 5.12 13.93 


y <- read.table("exchange.txt",header=F)
y <- y$V1
load("mcmc_exchange2")
load("stan_pd_vb.mf")
VA2apar <- read.table("VA2apar_exchange.txt", header=F)
VA2ameansd <- read.table("VA2ameansd_exchange.txt", header=F)
VA1diagbpar <- read.table("VA1diagbpar_exchange.txt", header=F)
VA1diagbmeansd <- read.table("VA1diagbmeansd_exchange.txt", header=F)


labels <- c(expression(alpha), expression(lambda), expression(psi))

n <- 945

pdf("C:\\Users\\statsll\\Dropbox\\Gaussian VA revised\\S&C_Manuscript\\images\\exchange.pdf", height=2.2, width=7)
par(mar=c(3.5,2.3,1,1))
par(mfrow=c(1,3))
par(mgp=c(2.7,0.7,0))
ymax <- c(3.7, 8.5, 1.1)
xmin <- c(-2.9, -1.5, 1.6)
xmax <- c(-0.8, 0.2, 7)
for (i in 1:3){
 temp <- density(mcmc[,i])
 plot(temp, main="", lwd=1.5, xlim=c(xmin[i],xmax[i]),
       ylab="", xlab=labels[i], ylim=c(0,ymax[i]), cex.lab=1.8, cex.axis=1.1)
 points(density(vb.mf[,i]),type="l",col="blue", lwd=1.5, lty=2)
 lines(temp$x, dnorm(temp$x,mean=VA1diagbmeansd[(n+i),1],sd=VA1diagbmeansd[(n+i),2]), 
      type='l', lty=1,col="blue", lwd=1.5, xlim=c(xmin[i],xmax[i]))
# lines(temp$x, dnorm(temp$x,mean=VA1bmeansd[(n+i),1],sd=VA1bmeansd[(n+i),2]), 
#      type='l', lty=1,col="green", lwd=1.5, xlim=c(xmin[i],xmax[i]))
 lines(temp$x, dnorm(temp$x,mean=VA2ameansd[(n+i),1],sd=VA2ameansd[(n+i),2]),
      type='l', lty=1,col="red", lwd=1.5, xlim=c(xmin[i],xmax[i]))
}
dev.off()



mcmch_mean <- apply(mcmc,2,mean)[4:(n+3)]
mcmch_sd <- apply(mcmc,2,sd)[4:(n+3)]


pdf("C:\\Users\\statsll\\Dropbox\\Gaussian VA revised\\S&C_Manuscript\\images\\exchange_h.pdf", height=3.7, width=8.5)
par(mar=c(3.5,4,1,1))
par(mgp=c(2.5,0.7,0))
#par(mfrow=c(2,1))
#plot(y,type="l")
plot(1:n, mcmch_mean, main="", lwd=1, type="l", ylim =c(-15,16),
       ylab=expression(b[t]), xlab="t", cex.lab=1.2, cex.axis=1)
lines(1:n, mcmch_mean-mcmch_sd, lwd=1, type="l", lty=2)
lines(1:n, mcmch_mean+mcmch_sd, lwd=1, type="l", lty=2)
lines(1:n, VA2ameansd[1:n,1], lwd=1, type="l", col="red")
lines(1:n, VA2ameansd[1:n,1] + VA2ameansd[1:n,2], lwd=1, type="l", lty=2, col="red")
lines(1:n, VA2ameansd[1:n,1] - VA2ameansd[1:n,2], lwd=1, type="l", lty=2, col="red")
dev.off()



par(mfrow=c(1,3))
for (i in 1:3){
plot(VA2apar[,i], type='l',col="red")
lines(VA1diagbpar[,i], type='l', col="blue")
lines(VA1bpar[,i], type='l', col="green")
}




##############################
# denmark exchange rate data #
##############################
require(Ecdat)
data(Garch)
ex <- Garch$dm
T <- length(ex)
tex <- log(ex[2:T]) - log(ex[1:(T-1)]) 
y <- 100*(tex -mean(tex))
plot(y,type="l")
n <- length(y)
# write.table(y, file="dm.txt",row.names=FALSE,col.names=FALSE)
data <- list(n=n, y=y)

# stan (mcmc) #
fit <- stan(file = 'model_svm2.stan', pars = c("alpha","lambda","psi","x"),
            data = data, iter = 50000, chains = 1, thin=5)
# print(fit, digits = 3)
la <- extract(fit, permuted = TRUE) # return a list of arrays 
mcmc <- cbind(la$alpha, la$lambda, la$psi, la$x)
(mcmc_mean <- apply(mcmc,2,mean))
(mcmc_sd <- apply(mcmc,2,sd))
# save(mcmc, file="mcmc_dm2")

par(mfrow=c(3,1))
par(mar=c(3,2,1,1))
for (i in 1:3){ plot(mcmc[,i],type='l') }

par(mfrow=c(1,3))
for (i in 1:3){ plot(density(mcmc[,i]))}

head(mcmc_mean)
head(mcmc_sd)

#  Elapsed Time: 2114.46 seconds (Warm-up)
#                2112.08 seconds (Sampling)
#                4226.54 seconds (Total)
> head(mcmc_mean)
[1] -1.6613624 -0.7689462  3.3776690 -6.0800261 -6.3059113 -6.4614228
> head(mcmc_sd)
[1] 0.1362314 0.1427573 0.3317822 2.8744907 2.7990974 2.7284060




# stan (compare mcmc, mean-field, fullrank)
model_svm <- stan_model(file='model_svm2.stan')

# mean-field #
system.time(fit.mf <- vb(model_svm, pars = c("alpha","lambda","psi","x"),
            data=data, algorithm="meanfield", output_samples=2000,seed=123))


   user  system elapsed 
  17.83    0.62   18.47 

print(fit.mf)
ex.mf <- extract(fit.mf, permuted = TRUE)
vb.mf <-  cbind(ex.mf$alpha, ex.mf$lambda, ex.mf$psi, ex.mf$x)
# save(vb.mf, file="stan_dm_vb.mf")
# Inference for Stan model: model_svm2.
1 chains, each with iter=2000; warmup=0; thin=1; 
post-warmup draws per chain=2000, total post-warmup draws=2000.

Inference for Stan model: model_svm2.
1 chains, each with iter=2000; warmup=0; thin=1; 
post-warmup draws per chain=2000, total post-warmup draws=2000.

         mean   sd  2.5%   25%   50%   75% 97.5%
alpha   -0.66 0.06 -0.78 -0.71 -0.66 -0.62 -0.55
lambda  -0.66 0.03 -0.71 -0.68 -0.66 -0.64 -0.61
psi     -4.64 1.06 -6.74 -5.39 -4.67 -3.91 -2.58


Approximate samples were drawn using VB(meanfield) at Fri Nov 11 13:30:54 2016.
We recommend genuine 'sampling' from the posterior distribution for final inferences!




# full-rank
system.time(fit.fr <- vb(model_svm, pars = c("alpha","lambda","psi","x"),
            data = data, algorithm="fullrank", output_samples=2000,seed=123))
print(fit.fr)
plot(fit.fr)
ex.fr <- extract(fit.fr, permuted = TRUE)
vb.fr <-  cbind(ex.mf$alpha, ex.mf$lambda, ex.mf$psi, ex.mf$x)

# This is Automatic Differentiation Variational Inference.

(EXPERIMENTAL ALGORITHM: expect frequent updates to the procedure.)

Gradient evaluation took 0 seconds
1000 iterations under these settings should take 0 seconds.
Adjust your expectations accordingly!

Begin eta adaptation.
Iteration:   1 / 250 [  0%]  (Adaptation)
Iteration:  50 / 250 [ 20%]  (Adaptation)
Iteration: 100 / 250 [ 40%]  (Adaptation)
Iteration: 150 / 250 [ 60%]  (Adaptation)
Iteration: 200 / 250 [ 80%]  (Adaptation)
Iteration: 250 / 250 [100%]  (Adaptation)
Error: stan::variational::advi::adapt_eta: All proposed step-sizes failed. Your model may be either severely ill-conditioned or misspecified.
Timing stopped at: 34.67 20.01 54.72 






VA2apar <- read.table("VA2apar_dm.txt", header=F)
VA2ameansd <- read.table("VA2ameansd_dm.txt", header=F)
VA1diagbpar <- read.table("VA1diagbpar_dm.txt", header=F)
VA1diagbmeansd <- read.table("VA1diagbmeansd_dm.txt", header=F)
load("mcmc_dm2")

labels <- c(expression(alpha), expression(lambda), expression(psi))
n <- dim(VA2ameansd)[1]-3

pdf("C:\\Users\\statsll\\Dropbox\\Gaussian VA revised\\S&C_Manuscript\\images\\dm.pdf", height=2.2, width=7)
par(mar=c(3.5,2.3,1,1))
par(mfrow=c(1,3))
par(mgp=c(2.7,0.7,0))
ymax <- c(5.6, 12.6, 1.8)
xmin <- c(-2.1, -1.2, 2.2)
xmax <- c(-1.2, -0.35, 4.8)
for (i in 1:3){
 temp <- density(mcmc[,i])
 plot(temp, main="", lwd=1.5, xlim=c(xmin[i],xmax[i]),
       ylab="", xlab=labels[i], ylim=c(0,ymax[i]), cex.lab=1.8, cex.axis=1.1)
 points(density(vb.mf[,i]),type="l",col="blue", lwd=1.5, lty=2)
 lines(temp$x, dnorm(temp$x,mean=VA1diagbmeansd[(n+i),1],sd=VA1diagbmeansd[(n+i),2]), 
      type='l', lty=1,col="blue", lwd=1.5, xlim=c(xmin[i],xmax[i]))
 lines(temp$x, dnorm(temp$x,mean=VA2ameansd[(n+i),1],sd=VA2ameansd[(n+i),2]),
      type='l', lty=1,col="red", lwd=1.5, xlim=c(xmin[i],xmax[i]))
}
dev.off()



mcmch_mean <- apply(mcmc,2,mean)[4:(n+3)]
mcmch_sd <- apply(mcmc,2,sd)[4:(n+3)]


pdf("C:\\Users\\statsll\\Dropbox\\Gaussian VA revised\\S&C_Manuscript\\images\\dm_h.pdf", height=3.6, width=10)
par(mar=c(3.5,3.5,1,1))
par(mgp=c(2.2,0.7,0))
plot(1:n, mcmch_mean, main="", lwd=1, type="l", ylim =c(-12,10),
       ylab=expression(b[t]), xlab="t", cex.lab=1.2, cex.axis=1)
lines(1:n, mcmch_mean-mcmch_sd, lwd=1, type="l", lty=2)
lines(1:n, mcmch_mean+mcmch_sd, lwd=1, type="l", lty=2)
lines(1:n, VA2ameansd[1:n,1], lwd=1, type="l", col="red")
lines(1:n, VA2ameansd[1:n,1] + VA2ameansd[1:n,2], lwd=1, type="l", lty=2, col="red")
lines(1:n, VA2ameansd[1:n,1] - VA2ameansd[1:n,2], lwd=1, type="l", lty=2, col="red")
dev.off()



