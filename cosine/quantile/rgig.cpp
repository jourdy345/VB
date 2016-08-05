// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <iostream>
#include <random>
#include <exception>
#include <stdexcept>
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

double _gig_mode(double lambda, double omega)
{
	if (lambda >= 1.)
		return (std::sqrt((lambda-1.)*(lambda-1.) + omega*omega)+(lambda-1.))/omega;
	else
		return omega / (std::sqrt((1.-lambda)*(1.-lambda) + omega*omega) + (1.-lambda));
}

void _rgig_ROU_noshift (double *res, int n, double lambda, double lambda_old, double omega, double alpha)
{
	double xm, nc;
	double ym, um;
	double s, t;
	double U,V, X;

	int i;
	int count = 0;

	t = 0.5 * (lambda-1.);
	s = 0.25 * omega;
	xm = _gig_mode(lambda, omega);

	nc = t*std::log(xm) - s*(xm+1./xm);

	ym = ((lambda+1.) + std::sqrt((lambda+1.)*(lambda+1.) + omega*omega))/omega;

	um = std::exp(0.5 *(lambda+1.)*std::log(ym)-s*(ym+1./ym)-nc);


	std::random_device generator;
	std::uniform_real_distribution<double> distribution(0.0, 1.0);
	for (i = 0; i < n; ++i) {
		do {
			++count;
			U = um * distribution(generator);
			V = distribution(generator);
			X = U/V;
		}
		while (std::log(V) > (t*std::log(X) - s*(X+1./X)-nc));

		res[i] = (lambda_old < 0.) ? (alpha / X) : (alpha * X);
	}

	return;
}

void _rgig_newapproach1(double *res, int n, double lambda, double lambda_old, double omega, double alpha)
{
	arma::vec A(3);
	double Atot;
	double k0;
	double k1, k2;

	double xm;
	double x0;
	double a;

	double U, V, X;
	double hx;

	int i;
	int count = 0;

	if (lambda >= 1. || omega > 1.)
		std::cout << "invalid parameters" << std::endl;
		std::exit(1);
		// throw std::invalid_argument("invalid parameters");

	xm = _gig_mode(lambda, omega);
	x0 = omega/(1.-lambda);

	k0 = std::exp((lambda-1.)*std::log(xm) - 0.5 * omega * (xm + 1./xm));
	A(0) = k0 * x0;

	if (x0 >= 2./omega) {
		k1 = 0.;
		A(1) = 0.;
		k2 = std::pow(x0, lambda-1.);
		A(2) = k2 * 2. * std::exp(-omega*x0/2.)/omega;
	}

	else {
		k1 = std::exp(-omega);
		A(1) = (lambda == 0.)
		  ? k1 * std::log(2./(omega*omega))
		  : k1 / lambda * (std::pow(2./omega, lambda) - std::pow(x0, lambda));

		  k2 = std::pow(2/omega, lambda-1.);
		  A(2) = k2 * 2 * std::exp(-1.)/omega;
	}

	Atot = arma::sum(A);
	std::random_device generator;
	std::uniform_real_distribution<double> distribution(0.0, 1.0);
	for (i = 0; i < n; ++i) {
		do {
			++count;
			V = Atot * distribution(generator);

			do {
				if (V <= A(0)) {
					X = x0 * V / A(0);
					hx = k0;
					break;
				}

				V -= A(0);
				if (V <= A(1)) {
					if (lambda == 0.) {
						X = omega * std::exp(std::exp(omega)*V);
						hx = k1 / X;
					}
					else {
						X = std::pow(std::pow(x0, lambda) + (lambda/k1 * V), 1./lambda);
						hx = k1 * std::pow(X, lambda - 1.);
					}
					break;
				}

				V -= A(1);
				a = (x0 > 2./omega) ? x0 : 2./omega;
				X = -2./omega * std::log(std::exp(-omega/2. * a) - omega/(2.*k2) * V);
				hx = k2 * std::exp(-omega/2. * X);
				break;
			} while(0);

			U = distribution(generator) * hx;

			if (std::log(U) <= (lambda-1.) * std::log(X) - omega/2. * (X+1./X)) {
				res[i] = (lambda_old < 0.) ? (alpha / X) : (alpha * X);
				break;
			}
		} while (1);
	}
	return;
}

void _rgig_ROU_shift_alt (double *res, int n, double lambda, double lambda_old, double omega, double alpha)
{
	double xm, nc;
	double s, t;
	double U, V, X;

	int i;
	int count = 0;
	double a, b, c;
	double p, q;
	double fi, fak;
	double y1, y2;

	double uplus, uminus;

	t = 0.5 * (lambda -1.);
	s = 0.25 * omega;

	xm = _gig_mode(lambda, omega);
	nc = t*std::log(xm)-s*(xm+1./xm);

	a = -(2. * (lambda + 1.) / omega + xm);
	b = (2. * (lambda - 1.) * xm / omega - 1.);
	c = xm;

	p = b - a*a/3.;
	q = (2. * a * a * a) / 27. - (a * b) / 3. + c;

	fi = std::acos(-q/(2.*std::sqrt(-(p*p*p)/27.)));
	fak = 2.*std::sqrt(-p/3.);
	y1 = fak * std::cos(fi/3.) - a/3.;
	y2 = fak * std::cos(fi/3. + 4./3.*M_PI) - a/3.;

	uplus = (y1 - xm) * std::exp(t*std::log(y1) - s*(y1 + 1./y1)-nc);
	uminus = (y2 - xm) * std::exp(t*std::log(y2) - s*(y2 + 1./y2) - nc);

	std::random_device generator;
	std::uniform_real_distribution<double> distribution(0.0, 1.0);
	for (i = 0; i < n; ++i) {
		do {
			++count;

			U = uminus + distribution(generator) * (uplus - uminus);
			V = distribution(generator);
			X = U/V + xm;
		}
		while ((X <= 0.) || ((std::log(V)) > (t * std::log(X) - s*(X + 1./X) - nc)));

		res[i] = (lambda_old < 0.) ? (alpha / X) : (alpha * X);
	}
	return;
}


arma::vec do_rgig(int n, double lambda, double chi, double psi)
/*-------------------------------------------------------------------*/
/* Draw sample from GIG distribution.                                */
/* 																																	 */
/* Parameters:																											 */
/*    n ....... sample size (positive integer)                       */
/*    lambda .. parameter for distribution                           */
/*    chi   ... parameter for distribition                           */
/*    psi   ... parameter for distribution                           */
/* Return:                                                           */
/*    random sample of size 'n'                                      */
/*-------------------------------------------------------------------*/
{
	const double ZTOL = (DOUBLE_EPS*10.0);
	double omega, alpha;
	/* check sample size */
	if (n <= 0) {
		std::cout << "sample size 'n' must be a positive integer" << std::endl;
		std::exit(1);
		// throw std::invalid_argument("sample size 'n' must be positive integer.");
	}

	/* check GIG parameters: */
	if ( !(R_FINITE(lambda) && R_FINITE(chi) && R_FINITE(psi)) ||
			 (chi < 0. || psi < 0)         ||
			 (chi == 0. && lambda <= 0.)   ||
			 (chi == 0. && lambda >= 0.)) {
		std::cout << "invalid parameters for GIG distribution" << std::endl;
		std::exit(1);
		// throw std::invalid_argument("invalid parameters for GIG distribution");
	}

	std::random_device generator;
	arma::vec res(n);
	if (chi < ZTOL) {
		/* special cases which are basically Gamma and Inverse Gamma distribution */
		if (lambda > 0.0) {
			for (int i = 0; i < n; ++i) {
				std::gamma_distribution<double> distribution(lambda, 2./psi);
				// res(i) =  rgamma(lambda, 2.0/psi);
				res(i) =  distribution(generator);
			} 
		}
		else {
			for (int i = 0; i < n; ++i) {
				std::gamma_distribution<double> distribution(-lambda, 2./psi);
				res(i) = 1.0/distribution(generator);
				// res(i) = 1.0/rgamma(-lambda, 2.0/psi);
			}
		}
	}

	else if (psi < ZTOL) {
		/* special cases which are basically Gamma and Inverse Gamma distribution */
		if (lambda > 0.0) {
			for (int i = 0; i < n; ++i) {
			  std::gamma_distribution<double> distribution(lambda, 2./chi);
			  res(i) = 1./distribution(generator);
			  // res(i) = 1.0/rgamma(lambda, 2.0/chi);
			}
		} else {
			for (int i = 0; i < n; ++i) {
				std::gamma_distribution<double> distribution(-lambda, 2./chi);
				res(i) = distribution(generator);
				// res(i) = rgamma(-lambda, 2.0/chi);
			}
		}
	}
	else {
		double lambda_old = lambda;
		if (lambda < 0.) lambda = -lambda;
		alpha = std::sqrt(chi / psi);
		omega = std::sqrt(psi*chi);

		do {
			if (lambda > 2. || omega > 3.) {
				/* Ratio of uniforms with shift by 'mode', alternative implementation */
				_rgig_ROU_shift_alt(res.memptr(), n, lambda, lambda_old, omega, alpha);
				break;
			}
			if (lambda >= 1.-2.25*omega*omega || omega > 0.2) {
				_rgig_ROU_noshift(res.memptr(), n, lambda, lambda_old, omega, alpha);
				break;
			}

			if (lambda >= 0. && omega > 0.) {
				_rgig_newapproach1(res.memptr(), n, lambda, lambda_old, omega, alpha);
				break;
			}

			std::cout << "parameters must satisfy lambda >= 0 and omega >= 0" << std::endl;
			std::exit(1);
			// throw invalid_argument("parameters must satisfy lambda>=0 and omega>=0");
		} while(0);
	}
	return res;
}

arma::vec rgig(int n, double lambda, double chi, double psi)
{
	arma::vec sexp_res = do_rgig(n, lambda, chi, psi);

	return sexp_res;
}

arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
   int ncols = sigma.n_cols;
   arma::mat Y = arma::randn(n, ncols);
   return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

double fRho(const double Yt, const double x, const double const1, const double tau, const arma::vec theta, const int J)
{
	double tmp = 0.;
	for (int i = 1; i < (J+1); ++i) {
		tmp += (std::exp(Yt * i) - std::exp(x * i)) * theta(i-1) * theta(i-1);
	}
	return std::exp(- tmp / (2. * tau) + const1 * (Yt - x) - x + Yt);
}

// [[Rcpp::export]]
Rcpp::List MCMCQuantile(arma::vec y, arma::vec x, arma::mat W, arma::mat varphi, int J, double p, int nSample, int burnIn, int thinIn, double A, double B, arma::vec muBeta, arma::mat SigmaBeta, double w0)
{
	int n = y.size();
	arma::mat SigmaBeta_inv = SigmaBeta.i();
	arma::colvec SigmaBeta_inv_muBeta = SigmaBeta_inv * muBeta;
	const double const1 = J * (J + 1.) / 4. - w0;
	const double nup    = (1. - 2. * p) / (p * (1. - p));
	const double taup2  = 2. / (p * (1. - p));
	std::random_device generator;

	/* initialize */
	arma::colvec beta(W.n_cols, arma::fill::randn);
	arma::colvec theta(J, arma::fill::randn);
	std::exponential_distribution<double> distribution(1.);
	// double gamma = distribution(generator);
	double gamma = 0.001;
	double tau   = 1./distribution(generator);
	const double const2 = 0.25 * taup2;
	arma::colvec u(n, arma::fill::zeros);
	std::uniform_real_distribution<double> unif_rand(0.0, 1.0);

	int count = 1;
	std::cout << "Burning in..." << std::endl;
	for (int t = 0; t < burnIn; ++t) {
		std::cout << "Burn-in count: " << count << std::endl;
		++count;
		/* update u */
		std::cout << ">" << std::endl;
		arma::vec tmp = y - W * beta - varphi * theta;
		tmp = (tmp % tmp) / taup2;
		std::cout << ">>" << std::endl;
		for (int i = 0; i < n; ++i) {
			u(i) = arma::as_scalar(rgig(1, 0.5, tmp(i), const2));
		}
		std::cout << ">>>" << std::endl;

		// u = rgig(n, 0.5, (tmp % tmp) / taup2, 0.25 *taup2);

		/* update beta */
		arma::mat W_halfu = arma::diagmat(1./arma::sqrt(u)) * W;
		arma::mat E_prev(J, J, arma::fill::zeros);
		std::cout << ">>>>" << count << std::endl;
		for (int k = 0; k < J; ++k) {
			E_prev(k, k) = std::exp((k+1) * gamma);
		}
		std::cout << ">>>>>" << count << std::endl;

		arma::mat E = E_prev / tau;
		arma::mat sigma_temp = (W_halfu.t() * W_halfu / taup2 + SigmaBeta_inv).i();
		arma::vec mu_temp    = sigma_temp * arma::vectorise(arma::sum((diagmat(((y - varphi * theta) - nup * u) / u) * W), 0) + SigmaBeta_inv_muBeta);
		beta = arma::vectorise(mvrnormArma(1, mu_temp, sigma_temp));
		std::cout << ">>>>>>" << count << std::endl;

		/* update theta */
		arma::mat varphi_halfu = arma::diagmat(1./arma::sqrt(u)) * varphi;
		sigma_temp = ((varphi_halfu.t() * varphi_halfu) / taup2 + E).i();
		mu_temp    = sigma_temp * arma::vectorise(arma::sum((arma::diagmat((y - W * beta - nup * u) / u) * varphi), 0));
		theta = arma::vectorise(mvrnormArma(1, mu_temp, sigma_temp));

		std::cout << ">>>>>>>" << count << std::endl;
		/* update tau */
		std::gamma_distribution<double> gam_dist(A + J/2., 1./(0.5 * arma::sum(arma::diagvec(E) % theta % theta) + B));
		tau = 1 / gam_dist(generator);

		std::cout << ">>>>>>>>" << count << std::endl;
		/* update gamma */
		double gamma_cand = distribution(generator);
		double rho = std::min(1., fRho(gamma_cand, gamma, const1, tau, theta, J));

		std::cout << ">>>>>>>>>" << count << std::endl;
		if (unif_rand(generator) < rho) {
			gamma = gamma_cand;
		}
		std::cout << ">>>>>>>>>>" << count << std::endl;
	}
	arma::mat uList(n, nSample+1);
	arma::mat betaList(W.n_cols, nSample);
	arma::mat thetaList(J, nSample);
	arma::vec tauList(nSample);
	arma::vec gammaList(nSample);
	uList.col(0) = u;
	betaList.col(0) = beta;
	thetaList.col(0) = theta;
	tauList(0) = tau;
	gammaList(0) = gamma;

	std::cout << "Thinning in..." << std::endl;
	count = 1;
	for (int j = 1; j < nSample; ++j) {
		std::cout << count << " samples collected..." << std::endl;
		++count;
		for (int t = 0; t < thinIn; ++t) {
			/* update u */
			arma::vec tmp = y - W * beta - varphi * theta;
			tmp = (tmp % tmp) / taup2;
			for (int i = 0; i < n; ++i) {
				u(i) = arma::as_scalar(rgig(1, 0.5, tmp(i), const2));
			}
			// u = rgig(n, 0.5, (tmp % tmp) / taup2, 0.25 *taup2);

			/* update beta */
			arma::mat W_halfu = arma::diagmat(1./arma::sqrt(u)) * W;
			arma::mat E_prev(J, J, arma::fill::zeros);
			for (int k = 0; k < J; ++k) {
				E_prev(k, k) = std::exp((k+1) * gamma);
			}
			arma::mat E = E_prev / tau;
			arma::mat sigma_temp = (W_halfu.t() * W_halfu / taup2 + SigmaBeta_inv).i();
			arma::vec mu_temp    = sigma_temp * arma::vectorise(arma::sum((diagmat(((y - varphi * theta) - nup * u) / u) * W), 0) + SigmaBeta_inv_muBeta);
			beta = arma::vectorise(mvrnormArma(1, mu_temp, sigma_temp));

			/* update theta */
			arma::mat varphi_halfu = arma::diagmat(1./arma::sqrt(u)) * varphi;
			sigma_temp = ((varphi_halfu.t() * varphi_halfu) / taup2 + E).i();
			mu_temp    = sigma_temp * arma::vectorise(arma::sum((arma::diagmat((y - W * beta - nup * u) / u) * varphi), 0));
			theta = arma::vectorise(mvrnormArma(1, mu_temp, sigma_temp));

			/* update tau */
			std::gamma_distribution<double> gam_dist(A + J/2., 1./(0.5 * arma::sum(arma::diagvec(E) % theta % theta) + B));
			tau = 1 / gam_dist(generator);

			/* update gamma */
			double gamma_cand = distribution(generator);
			double rho = std::min(1., fRho(gamma_cand, gamma, const1, tau, theta, J));

			if (unif_rand(generator) < rho) {
				gamma = gamma_cand;
			}

			// /* update u */
			// arma::vec tmp = y - W * betaList.col(t-1) - varphi * thetaList.col(t-1);
			// tmp = (tmp % tmp) / taup2;

			// for (int i = 0; i < n; ++i) {
			// 	uList(i, t) = arma::as_scalar(rgig(1, 0.5, tmp(i), const2));
			// }
			// // uList.col(t) = rgig(n, 0.5, (tmp % tmp) / taup2, 0.25 *taup2);

			// /* update beta */
			// arma::mat W_halfu = arma::diagmat(arma::sqrt(uList.col(t-1))) * W;
			// arma::mat E_prev(J, J, arma::fill::zeros);
			// for (int k = 0; k < J; ++k) {
			// 	E_prev(k, k) = (k+1) * gammaList(t-1);
			// }
			// arma::mat E = E_prev / tauList(t-1);
			// arma::mat sigma_temp = (W_halfu.t() * W_halfu + SigmaBeta_inv).i();
			// arma::vec mu_temp = sigma_temp * arma::vectorise(arma::sum((diagmat(((y - varphi * thetaList.col(t-1)) - nup * uList.col(t-1)) / uList.col(t-1)) * W), 0) + SigmaBeta_inv_muBeta);
			// betaList.col(t) = arma::vectorise(mvrnormArma(1, mu_temp, sigma_temp));

			// /* update theta */
			// arma::mat varphi_halfu = arma::diagmat(arma::sqrt(uList.col(t-1))) * varphi;
			// sigma_temp = ((varphi_halfu.t() * varphi_halfu) / taup2 + E).i();
			// mu_temp = sigma_temp * arma::vectorise(arma::sum((arma::diagmat((y - W * betaList.col(t-1) - nup * uList.col(t-1)) / uList.col(t-1)) * varphi), 0));
			// thetaList.col(t) = arma::vectorise(mvrnormArma(1, mu_temp, sigma_temp));

			// /* update tau */
			// std::gamma_distribution<double> gam_dist(A + J/2., 1./(0.5 * arma::sum(arma::diagvec(E) % thetaList.col(t-1) % thetaList.col(t-1)) + B));
			// tauList(t) = 1 / gam_dist(generator);
			
			// double gamma_cand = distribution(generator);
			// double rho = fRho(gamma_cand, gamma, const1, tauList(t-1), thetaList.col(t-1), J);

			// if (unif_rand(generator) < rho) {
			// 	gammaList(t) = gamma_cand;
			// } else {
			// 	gammaList(t) = gammaList(t-1);
			// }
		}
		uList.col(j) = u;
		thetaList.col(j) = theta;
		betaList.col(j) = beta;
		tauList(j) = tau;
		gammaList(j) = gamma;
	}
	return Rcpp::List::create(Rcpp::Named("u") = uList,
										        Rcpp::Named("theta") = thetaList,
										        Rcpp::Named("beta") = betaList,
										        Rcpp::Named("tau") = tauList,
										        Rcpp::Named("gamma") = gammaList);
}	











