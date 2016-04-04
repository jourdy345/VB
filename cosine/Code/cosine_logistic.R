lambda_xi = function(x) -tanh(x/2)/(4*x)
Psi_xi = function(x) x/2 - log(1 + exp(x)) + x*tanh(x/2)/4
tr = function(X) sum(diag(X))
Qj = function(mu, sigma2, j) {
  t1 = 0.5 * sigma2 * j^2
  t2 = mu * j
  t3 = mu/sqrt(sigma2)
  t4 = sqrt(sigma2) * j
  (exp(term1 + term2) * (1 - pnorm(-(term3 + term4)))) + (exp(term1 - term2) * (1 - pnorm(term3 - term4)))
}
S1 = function(mu, sigma2, w0) -w0 * ((sqrt(sigma2 * 2 / pi) * exp(-mu^2 / (2 * sigma2))) + (mu * (1 - 2 * pnorm(-mu/sqrt(sigma2)))))
S2 = function(mupsi2, sigpsi2, w0, J, rtq, stq, sigtq, mutq) (-0.5 * rtq/stq * sum((diag(sigtq) + mutq^2) * Qj(mupsi, sigpsi2, 1:J))) - (J*(J+1)/(4*w0) * S1(mupsi, sigpsi2, w0))
LB = function(mutq, varphi, xi, sigtq, rtq, stq, rt0, st0, w0, mupsi, sigpsi2) {
  J = dim(varphi)[2]
  t1 = as.vector(varphi %*% mutq)
  lambda = lambda_xi(xi)
  half_rtq = rtq / 2
  half_stq = stq / 2
  half_rt0 = rt0 / 2
  half_st0 = st0 / 2
  ratioq = rtq / stq
  # varphi = sqrt(2)*cos(outer(x,pi*(1:J)))
  crossprod(t1, y - 0.5) + sum(lambda * diag(varphi %*% tcrossprod(sigtq, varphi))) + sum(t1^2 * lambda) + sum(Psi_xi(xi))
  + 0.5 * J * (digamma(half_rtq) - log(half_stq) - log(2*pi)) + S1(mupsi, sigpsi2, w0) + S2(mupsiq, sigpsi2, w0, J, rtq, stq, sigtq, mutq)
  + half_rt0 * log(half_st0) - lgamma(half_rt0) + half_rt0 * (digamma(half_rtq) - log(half_stq)) - half_st0 * ratioq - log(w0/2)
  + half_rtq + lgamma(half_rtq) - half_rtq * digamma(half_rtq) + 0.5 * (J * (1 + log(2*pi)) + determinant(sigtq)$modulus[1] + log(2*pi*sigpsi2) + 1)
}


VB = function(y, x, rt0, st0, w0, J, tol = 1.0e-6) {
  varphi = sqrt(2)*cos(outer(x,pi*(1:J)))

  # Initialize variational parameters
  rtq =
  stq = 
  ratioq = rtq/stq
  sigtq = 
  mutq = 
  xi = 
  sigpsi2 =
  mupsi = 

  while(dif > tol) {
    # Update theta
    sigtq = solve(crossprod(varphi, varphi * lambda_xi(xi)) - 0.5 * ratioq * diag(Qj(mupsi, sigpsi2, 1:J)))
    mutq = sigtq %*% crossprod(varphi, -2 * y - 1)

    # Update tau
    stq = st0 + sum((diag(sigtq) + mutq^2) * Qj(mupsi, sigpsi2, 1:J))

    # Update psi
    sigpsi2 = 
  }

}