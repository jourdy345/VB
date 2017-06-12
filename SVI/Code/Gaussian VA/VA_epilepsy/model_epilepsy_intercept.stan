data {
  int<lower=0> n; // number of subjects 
  int<lower=0> N; // number of obs 
  int<lower=0> k; // number of fixed effects
  int<lower=0> y[N]; // responses
  matrix[N,k] X; // fixed effects covariates
  vector[N] Z; // random effects covariates 
  int<lower=1> startindex[n];
  int<lower=1> endindex[n];
}
parameters {
  vector[k] beta;
  real zeta;
  vector[n] b;
}
model {
  vector[N] prob;
  zeta ~ normal(0, 10);  // vectorized form, 10 is the standard deviation
  beta ~ normal(0, 10);  // vectorized form  
  for (i in 1:n) {
   b[i] ~ normal(0, exp(zeta)); 
   for (j in startindex[i]:endindex[i]){
    prob[j] <- dot_product(X[j,], beta) + Z[j]* b[i];
   }
  }
  y ~ poisson_log(prob);
}
