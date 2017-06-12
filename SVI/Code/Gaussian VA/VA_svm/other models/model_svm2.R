model {

  for (t in 1:n){
   y[t] ~ dnorm(0, py[t])
   py[t] <- exp(lambda)*exp(alpha*x[t])
  }

  x[1] ~ dnorm(0,p1)
  p1 <- 1- pow(phi,2)

  for (t in 2:n){
   x[t] ~ dnorm(mu[t-1], 1)
   mu[t-1] <- phi*x[t-1]
  }

  ph1 <- (exp(psi) - 1)/(exp(psi) + 1)

  psi ~ dnorm(0, 0.1)
  lambda ~ dnorm(0, 0.1)
  alpha ~ dnorm(0, 0.1)
}
