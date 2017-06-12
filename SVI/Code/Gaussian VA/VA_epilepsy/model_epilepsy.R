model{

for (i in 1:n) {
 b[i,1:p] ~ dmnorm(mu[], P[,])
 for (j in startindex[i]:endindex[i]){
  y[j] ~ dpois(m[j])
  log(m[j]) <- inprod(X[j,], beta[]) + inprod(Z[j,], b[i,])
 }
}

for (l in 1:k){
 beta[l] ~ dnorm(0, 0.01)
}

# P: precision matrix, Sigma: covariance matrix
P[1:p,1:p] <- inverse(Sigma[,])

# Sigma = WW^T
for (g in 1:p){
 mu[g] <- 0
 for (h in 1:p){
  Sigma[g,h] <- inprod(W[g,], W[h,])
 }
}
 
W[1,1] <- exp(zeta[1])
W[1,2] <- 0
W[2,1] <- zeta[2]
W[2,2] <- exp(zeta[3])

for (s in 1:3){
 zeta[s] ~ dnorm(0, 0.01)
}

} 
