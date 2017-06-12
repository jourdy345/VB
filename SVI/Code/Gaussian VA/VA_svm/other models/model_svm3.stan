data {
  int<lower=0> n; // number of subjects 
  vector[n] y; // responses
}
parameters {
  vector[n] h;
  real alpha;
  real lambda;
  real psi;
}
transformed parameters{
  real phi;
  phi <- (exp(psi)-1)/(exp(psi)+1); 
} 
model {
  alpha ~ normal(0, sqrt(10) );
  lambda ~ normal(0, sqrt(10) );
  psi ~ normal(0, sqrt(10) );

  h[1] ~ normal(lambda, exp(alpha/2)/sqrt(1-phi^2) );

  for (t in 2:n){
   h[t] ~ normal(phi*h[t-1] + (1-phi)*lambda, exp(alpha/2) );
  }

  for (t in 1:n){
   y[t] ~ normal(0, exp(h[t]/2));
  }

}
