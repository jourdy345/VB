# Package for Dirichlet distribution
## Usage: rdirichlet, ddirichlet
library(gtools)

# Acceptance-Rejection algorithm for q(m_j)
m_sampling <- function(n, Cm, nm1, nm0) {
  f_div_g <- function(x) exp(-Cm * x) * x^(nm1-1) * (1-x)^(nm0-1) / dbeta(x, nm1, nm0)
  M <- optimize(f_div_g, interval = c(0, 1), maximum = TRUE)$objective
  temp <- c()
  while (length(temp) < n) {
    u <- runif(1)
    x <- rbeta(1, nm1, nm0)
    if (M * u < f_div_g(x)) temp <- c(temp, x)
  }
  temp
}

# Normalizing constant for q(m_j)
I_mj <- function(a, alpha, beta) {
  x <- rbeta(10000, alpha, beta)
  mean(exp(-a*x)) * beta(alpha, beta)
}

# Expectation of m_j with respect to posterior q(m_j)
Emj <- function(Cm, nm1, nm0) {
  I_mj(Cm, nm1 + 1, nm0) / I_mj(Cm, nm1, nm0)
}

# Vectorise Emj
Emj_v <- function(Cm, nm1, nm0) {
  m <- length(Cm)
  temp <- c()
  for (j in 1:m) {
    temp <- c(temp, Emj(Cm[j], nm1, nm0))
  }
  temp
}

# Expectation of log(Gamma(s_j)) with respect to posterior q(s_j)
E_lgamma_s <- function(a_s, beta_s) {
  mean(lgamma(rgamma(10000, shape = a_s, scale = beta_s)))
}

# Expectation of log(Gamma(s_j * m_j)) with respect to posterior q(s_j) & q(m_j)
E_lgamma_sm <- function(a_s, beta_s, Cm, nm1, nm0) {
  mean(lgamma(rgamma(10000, shape = a_s, scale = beta_s) * m_sampling(10000, Cm, nm1, nm0)))
}

# Expectation of log(Gamma(s_j - s_j * m_j - 1)) with respect to posterior q(s_j) & q(m_j)
E_lgamma_s_sm_1 <- function(a_s, beta_s, Cm, nm1, nm0) {
  s <- rgamma(10000, shape = a_s, scale = beta_s)
  m <- m_sampling(10000, Cm, nm1, nm0)
  mean(lgamma(s - s * m - 1))
}

# Lower bound
LB <- function(y, phi, a, a_s, beta_s, Cm, nm1, nm0, N, M) {
  # temp1
  temp1 <- 0
  for (i in 1:N) {
    for (j in 1:M) {
      temp1 <- temp1 + phi[i, j] * E_lgamma_s(a_s, beta_s[j]) - E_lgamma_sm(a_s, beta_s[j], Cm[j], nm1, nm0) - E_lgamma_s_sm_1(a_s, beta_s[j], Cm[j], nm1, nm0) + (a_s * beta_s[j] * Emj(Cm[j], nm1, nm0) - 1) * log(y[i]) + (a_s * beta_s[j] * (1 - Emj(Cm[j], nm1, nm0)) - 1) * log(1 - y[i])
    }
  }
  
  # temp2
  temp2 <- 0
  for (j in 1:M) {
    for (i in 1:N) {
      temp2 <- temp2 + digamma(phi[i, j] + a) - digamma(sum(phi) + M * a)
    }
  }

  # temp3
  temp3 <- 0
  for (j in 1:M) {
    temp3 <- temp3 - a_s * beta_s[j] + Cm[j] * Emj(Cm[j], nm1, nm0) + log(I_mj(Cm[j], nm1, nm0))
  }

  # temp4
  temp4 <- 0
  for (j in 1:M) {
    temp4 <- temp4 - a_s * beta_s[j] - 1 / (a_s - 1)
  }
  temp4 <- temp4 - M * lgamma(a_s) + sum(phi * log(phi))

  # Final lower bound
  temp1 - temp2 + temp3 - temp4
}

# Variational approximation for Beta Mixture
## y: N x 1 vector
## a, a_s, b_s, nm1, nm0: scalar
## M: number of mixtures (scalar)
BMVB <- function(y, M, a, a_s, b_s, nm1, nm0) {
  # Initialize variational parameters
  N <- length(y)
  phi <- rdirichlet(N, rep(a, M))
  beta_s <- rep(b_s, M)
  LB_old <- -Inf

  # Coordinate ascent
  while (TRUE) {
    ## Cm
    Cm <- -apply(phi * log(y / (1 - y)), 2, sum) * (a_s * beta_s)
    print(paste0('Cm ', Cm))
    ## beta_s
    beta_s <- 1 / (-((Emj_v(Cm, nm1, nm0) * apply(phi * log(y), 2, sum)) + ((1 - Emj_v(Cm, nm1, nm0)) * apply(phi * log(1 - y), 2, sum))) + 1/b_s)
    print(paste0('beta_s ', beta_s))
    ## phi_ij
    for (j in 1:M) {
      for (i in 1:N) {
        phi[i, j] <- (y[i] / (1 - y[i]))^(a_s * beta_s[j] * Emj(Cm[j], nm1, nm0)) * (1 - y[i])^(a_s * beta_s[j]) * ((sum(phi[, j]) + a) / (sum(phi) + M * a)) / (y[i] * (1 - y[i]))
      }
    }

    # Monitor convergence
    LB_new <- LB(y, phi, a, a_s, beta_s, Cm, nm1, nm0, N, M)
    if (abs(LB_new - LB_old) < 0.000001) break
    else LB_old <- LB_new
  }
  list('phi' = phi, 'beta_s' = beta_s, 'Cm' = Cm, 'LB' = LB_new)
}

# Simulation
M <- 3
N <- 300
## Hyperpriors
a <- 3
a_s <- 3
b_s <- 10^2
nm1 <- 2
nm0 <- 2
## Sampling priors from hyperpriors
s <- rgamma(M, shape = a_s, scale = b_s)
m <- rbeta(M, nm1, nm0)
lambda <- rdirichlet(1, rep(a, M))
## Generating mixture data (Avoid 'for' loops for performance reasons)
components <- sample(1:M, prob = lambda, size = N, replace = TRUE)
y <- rbeta(n = N, shape1 = s[components] * m[components], shape2 = s[components] * (1 - m[components]))
# Beta mixture Variational Bayes
BMVB(y, M, a, a_s, b_s, nm1, nm0)
