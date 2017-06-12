data {
  int<lower=0> n; // number of subjects 
  vector[n] y; // responses
}
parameters {
  vector[n] x;
  real psi;
  real lambda;
  real alpha;
}
transformed parameters{
  real phi;
  phi <- 1/(exp(-psi)+1); 
} 
model {
  alpha ~ normal(0, sqrt(10) );
  psi ~ normal(0, sqrt(10) );
  lambda ~ normal(0, sqrt(10) );

  x[1] ~ normal(0, 1/sqrt(1-phi^2) );

  for (t in 2:n){
   x[t] ~ normal(phi*x[t-1], 1);
  }

  for (t in 1:n){
   y[t] ~ normal( 0, exp(lambda/2 + exp(alpha)*x[t]/2) );
  }

}

