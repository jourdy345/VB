library(Rcpp)
library(RcppArmadillo)
library(quantmod)
sourceCpp('rgig.cpp')


nasdaq <- new.env()
status <- getSymbols('GOOGL', env = nasdaq, from = as.Date('2005-01-01'))
status <- getSymbols('AAPL', env = nasdaq, from = as.Date('2005-01-01'))

apple  <- get('AAPL', envir = nasdaq)
google <- get('GOOGL', envir = nasdaq)

x <- as.vector(apple[,4])
y <- as.vector(google[,4])

wrapMCMC <- function(y, x, W, prior, J, p, nSample, burnIn, thinIn) {
  if (!is.matrix(W)) W <- as.matrix(W)
  A <- prior$A
  B <- prior$B
  muBeta <- prior$muBeta
  SigmaBeta <- prior$SigmaBeta
  w0 <- prior$w0

  if (any(x < 0 | x > 1)) {
    xmax <- max(x)
    xmin <- min(x)
    x_adj <- (x - xmin) / (xmax - xmin)
    varphi <- sqrt(2/(xmax - xmin))*cos(outer(x_adj, pi*(1:J)))
  } else {
    varphi <- sqrt(2) * cos(outer(x, pi(1:J)))
  }

  # arma::vec y, arma::vec x, arma::mat W, arma::mat varphi, int J, double p, int nSample, int burnIn, int thinIn, double A, double B, arma::vec muBeta, arma::mat SigmaBeta, double w0
  res <- MCMCQuantile(y, x, W, varphi, J, p, nSample, burnIn, thinIn, A, B, muBeta, SigmaBeta, w0)
  res
}

n <- length(x)
W <- matrix(1, nr = n, nc = 1)
A <- 0.01
B <- 0.01
muBeta <- rep(0, ncol(W))
SigmaBeta <- diag(1, ncol(W))
w0 <- 1
prior <- list(A = A, B = B, muBeta = muBeta, SigmaBeta = SigmaBeta, w0 = w0)
J <- 20
p <- 0.3
fit <- wrapMCMC(y, x, W, prior, J, p, 4000, 4000, 5)