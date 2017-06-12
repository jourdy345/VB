model {

for (t in 1:n){
 y[t] ~ dnorm(0, py[t])
 py[t] <- exp(-eps)*exp(-h[t])
}

h[1] ~ dnorm(0,ph1)
ph1 <- (1- pow(tphi,2))*exp(-v)
for (t in 2:n){
 h[t] ~ dnorm(muh[t-1], ph)
 muh[t-1] <- tphi*h[t-1]
}

ph <- exp(-v)
tphi <- exp(phi)/(1+exp(phi))
phi ~ dnorm(0, 1)
v ~ dnorm(0, 1)
eps ~ dnorm(0, 1)
}
