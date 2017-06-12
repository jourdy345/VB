model{

for (i in 1:n) {
 b[i] ~ dnorm(0, P)
 for (j in startindex[i]:endindex[i]){
  y[j] ~ dbern(prob[j])
  logit(prob[j]) <- inprod(X[j,], beta[]) + Z[j]*b[i]
 }
}

for (l in 1:k){
 beta[l] ~ dnorm(0, 0.01)
}

P <- pow(W,-2)
W <- exp(zeta)
zeta ~ dnorm(0, 0.01)

} 
