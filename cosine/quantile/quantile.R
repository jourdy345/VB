# f <- function(x) exp(2*x)*exp(-4*exp(2*x))
# g <- function(x) exp(-x)


# burnin <- 10000
# X <- c(2)
# for (t in 2:burnin) {
# 	y_cand <- rexp(1)
# 	rho <- f(y_cand) * g(X[t-1]) / (f(X[t-1]) * g(y_cand))
# 	X[t] <- X[t-1] + (y_cand - X[t-1]) * (runif(1) < rho)
# }

# thinin <- 40
# N <- 8000
# X <- X[length(X)]
# sample <- c()
# for (j in 1:N) {
# 	for(i in 2:thinin) {
# 		y_cand <- rexp(1)
# 		rho <- f(y_cand) * g(X[i-1]) / (f(X[i-1]) * g(y_cand))
# 		X[i] <- X[i-1] + (y_cand - X[i-1]) * (runif(1) < rho)
# 	}
# 	sample <- c(sample, X[length(X)])
# }

require('GIGrvg')
require('mvtnorm')
fRho <- function(Yt, x, const1, tau, theta, J) {
	exp(const1 * (Yt - x) - sum((exp(Yt*(1:J)) - exp(x*(1:J))) * theta^2) / (2 * tau) - x + Yt)
}
MCMC <- function(y, x, W, prior, J, p, n_sample, burnin, thinin) {
	n <- length(y)
	if (!is.matrix(W)) W <- as.matrix(W)
	A <- prior$A
	B <- prior$B
	muBeta <- prior$muBeta
	SigmaBeta <- prior$SigmaBeta
	w0 <- prior$w0
	SigmaBeta_inv <- solve(SigmaBeta)
	SigmaBeta_inv_muBeta <- as.vector(SigmaBeta_inv %*% muBeta)
	xmax <- max(x)
	xmin <- min(x)
	x_adj <- (x - xmin) / (xmax - xmin)
	varphi <- sqrt(2/(xmax - xmin))*cos(outer(x_adj, pi*(1:J)))

	const1 <- J * (J + 1) / 4 - w0
	#initialize
	beta <- rnorm(ncol(W))
	theta <- rnorm(J)
	gamma <- rexp(1)
	tau   <- 1/rexp(1)

	nup <- (1-2*p)/(p*(1-p))
	taup2 <- 2/(p*(1-p))



	for (t in 1:burnin) {
		# update u
			u <- rgig(n, 0.5, (y - as.vector(W %*% beta) - as.vector(varphi %*% theta))^2 / taup2, 0.25 * taup2)

			# update beta
			W_halfu <- diag(sqrt(u)) %*% W
			E_prev <- diag((1:J)*gamma)
			E      <- E_prev / tau
			sigma_temp <- solve(crossprod(W_halfu)/taup2 + SigmaBeta_inv)
			mu_temp    <- sigma_temp %*% (colSums(diag((y - as.vector(varphi %*% theta) - nup * u) / u) %*% W) + SigmaBeta_inv_muBeta)
			beta <- as.vector(rmvnorm(1, mean = mu_temp, sigma = sigma_temp))

			# update theta
			varphi_halfu <- diag(sqrt(u)) %*% varphi
			sigma_temp <- solve(crossprod(varphi_halfu)/taup2 + E)
			mu_temp    <- sigma_temp %*% colSums(diag((y - as.vector(W %*% beta) - nup * u) / u) %*% varphi)
			theta <- as.vector(rmvnorm(1, mean = mu_temp, sigma = sigma_temp))

			# update tau
			tau <- 1/rgamma(1, A+J/2, rate = 0.5 * sum(diag(E) * theta^2) + B)

			gamma_cand <- rexp(1)
			# rho <- exp(const1 * gamma_cand - sum(exp((1:J)*gamma_cand) * theta^2) / (2 * tau)) * dexp(gamma) / (exp(const1 * gamma - sum(exp((1:J)*gamma) * theta^2) / (2 * tau)) * dexp(gamma_cand))
			rho <- fRho(gamma_cand, gamma, const1, tau, theta, J)

			gamma <- gamma + (gamma_cand - gamma) * (runif(1) < rho)

	}

	u <- list(u)
	theta <- list(theta)


	for (j in 1:n_sample) {
		for (t in 2:thinin) {
			# update u
			u[[t]] <- rgig(n, 0.5, (y - as.vector(W %*% beta[[t-1]]) - as.vector(varphi %*% theta[[t-1]]))^2 / taup2, 0.25 * taup2)

			# update beta
			W_halfu <- diag(sqrt(u[[t-1]])) %*% W
			E_prev <- diag((1:J)*gamma[t-1])
			E      <- E_prev / tau[t-1]
			sigma_temp <- solve(crossprod(W_halfu)/taup2 + SigmaBeta_inv)
			mu_temp    <- sigma_temp %*% (colSums(diag((y - as.vector(varphi %*% theta[[t-1]]) - nup * u[[t-1]]) / u[[t-1]]) %*% W) + SigmaBeta_inv_muBeta)
			beta[[t]] <- as.vector(rmvnorm(1, mean = mu_temp, sigma = sigma_temp))

			# update theta
			varphi_halfu <- diag(sqrt(u[[t-1]])) %*% varphi
			sigma_temp <- solve(crossprod(varphi_halfu)/taup2 + E)
			mu_temp    <- sigma_temp %*% colSums(diag((y - as.vector(W %*% beta[[t-1]]) - nup * u[[t-1]]) / u[[t-1]]) %*% varphi)
			theta[[t]] <- as.vector(rmvnorm(1, mean = mu_temp, sigma = sigma_temp))

			# update tau
			tau[t] <- 1/rgamma(1, A+J/2, rate = 0.5 * sum(diag(E) * theta[[t-1]]^2) + B)

			gamma_cand <- rexp(1)
			# rho <- exp(const1 * gamma_cand - sum(exp((1:J)*gamma_cand) * theta[[t-1]]^2) / (2 * tau[t-1])) * dexp(gamma[t-1]) / (exp(const1 * gamma[t-1] - sum(exp((1:J)*gamma[t-1]) * theta[[t-1]]^2) / (2 * tau[t-1])) * dexp(gamma_cand))
			rho <- fRho(gamma_cand, gamma[t-1], const1, tau[t-1], theta[[t-1]], J)
			gamma[t] <- gamma[t-1] + (gamma_cand - gamma[t-1]) * (runif(1) < rho)
		}
	}
	u <- do.call('cbind', u)
	theta <- do.call('cbind', theta)
	list(u = u, theta = theta, gamma = gamma, tau = tau)
}


