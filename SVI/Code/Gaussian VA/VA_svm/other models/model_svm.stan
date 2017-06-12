data {
  int<lower=0> n; // number of subjects 
  vector[n] y; // responses
}
parameters {
  vector[n] h;
  real<lower=0,upper=1> phistar;
  real mu;
  real<lower=0> v;
}
transformed parameters{
  real phi;
  phi <- 2*phistar-1; 
} 
model {
  v ~ inv_gamma(2.5, 0.025);
  phistar ~ beta(20,1.5);
  mu ~ normal(0,sqrt(10));

  h[1] ~ normal(mu, sqrt(v/(1-phi^2)));

  for (t in 2:n){
   h[t] ~ normal(phi*h[t-1] + (1-phi)*mu, sqrt(v));
  }

  for (t in 1:n){
   y[t] ~ normal(0, exp(h[t]/2));
  }

}
