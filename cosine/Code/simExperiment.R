# require(Rcpp)
# require(RcppArmadillo)
# sourceCpp('auxiliary.cpp')
# simExperiment = function(FUN, J, draw = TRUE, intercept = TRUE) {

#   dnormOverPnorm <- function(x) sqrt(2 / pi) * exp(-x^2 / 2) / (erf(x / sqrt(2)) + 1)
#   dnormOverPnorm_Laurent <- function(x) -x - 1 / x + 2 / x^3 - 10 / x^5
#   dnormOverPnormMinusOne <- function(x) exp(-x^2 / 2) / (sqrt(pi / 2) * (erf(x / sqrt(2)) + 1) - 1)

#   tr = function(X) sum(diag(X))
#   E_foldedNormal = function(mu, sigma2, sigma) {
#     (sigma * sqrt(2 / pi) * exp(-mu^2/(2 * sigma2))) + (mu * (1 - 2 * pnorm(-mu / sigma)))
#   }
#   MGF_foldedNormal = function(sigma2, sigma, mu, t) {
#     (exp( (0.5 * sigma2 * t^2) + (mu * t)) * (1 - pnorm(-mu/sigma - (sigma * t)))) + (exp((0.5 * sigma2 * t^2) - (mu * t)) * (1 - pnorm(mu/sigma - (sigma * t))))
#   }
#   der_Qj_mu = function(sigma2, sigma, mu, j) {
#     (exp((0.5 * sigma2 * j^2) + (mu * j)) * dnorm(-mu/sigma - sigma * j) / sigma) + (j * exp((0.5 * sigma2 * j^2) + (mu * j)) * (1 - pnorm((-mu/sigma) - (sigma * j)))) - (exp(0.5 * sigma2 * j^2 - mu * j) * dnorm(mu/sigma - sigma * j) / sigma) - (j * exp((0.5 * sigma2 * j^2) - (mu * j)) * (1 - pnorm(mu/sigma - sigma * j)))
#   }
#   der_Qj_sigma2 = function(sigma2, sigma, mu, j) {
#     (-exp(0.5 * sigma2 * j^2 + mu * j) * dnorm(-mu/sigma - sigma * j) * ((mu/(2 * sigma^3)) - (j / (2 * sigma)))) + (0.5 * j^2 * exp(0.5 * sigma2 * j^2 + mu * j) * (1 - pnorm(-mu/sigma - (sigma * j)))) + (exp((0.5 * sigma2 * j^2) - (mu * j)) * dnorm(mu/sigma - (sigma * j)) * ((mu/(2 * sigma^3)) + (j/(2 * sigma)))) + (0.5 * j^2 * exp((0.5 * sigma2 * j^2) - (mu * j)) * (1 - pnorm(mu/sigma - sigma * j)))
#   }

#   der_S1_mu = function(sigma2, sigma, mu, w0) {
#     -w0 * ( (-mu/sigma * sqrt(2/pi) * exp(-mu^2 / (2 * sigma2))) + (1 - 2 * pnorm(-mu/sigma)) + (2*mu/sigma * dnorm(-mu/sigma)) )
#   }

#   der_S1_sigma2 = function(sigma2, sigma, mu, w0) {
#     -w0 * (((1/(2 * sigma)) + (0.5 * mu^2 / sigma^3) * sqrt(2/pi) * exp(-mu^2 / (2 * sigma2))) - (mu^2 / sigma^3 * dnorm(-mu/sigma)))
#   }

#   der_S2_mu = function(sigma2, sigma, sigtq, mutq, mu, rsoverssq, rtoverstq, w0, J) {
#     ((-0.25 * J * (J + 1) / w0) * der_S1_mu(sigma2, sigma, mu, w0)) - (0.5 * rsoverssq * rtoverstq * sum(((diag(sigtq) + mutq^2) * der_Qj_mu(sigma2, sigma, mu, (1:J)))))
#   }

#   der_S2_sigma2 = function(sigma2, sigma, sigtq, mutq, mu, rsoverssq, rtoverstq, w0, J) {
#     ((-0.25 * J * (J + 1) / w0) * der_S1_sigma2(sigma2, sigma, mu, w0)) - (0.5 * rsoverssq * rtoverstq * sum(((diag(sigtq) + mutq^2) * der_Qj_sigma2(sigma2, sigma, mu, (1:J)))))
#   }


#   LB = function(y, x, W, sigbq, varphi, sigtq, mubq, mutq, sigpsiq2, mupsiq, rsq, ssq, rtq, stq, rs0, ss0, rt0, st0, sigb0, mub0, w0) {
#     J = length(mutq)
#     p = length(mubq)
#     WtW = crossprod(W)
#     varphitvarphi = crossprod(varphi)
#     tmp = W %*% mubq + varphi %*% mutq
#     sigpsiq = sqrt(sigpsiq2)
#     rt0_half = 0.5 * rt0
#     st0_half = 0.5 * st0
#     rtq_half = 0.5 * rtq
#     stq_half = 0.5 * stq
#     rtoverstq = rtq / stq
#     rs0_half = 0.5 * rs0
#     ss0_half = 0.5 * ss0
#     rsq_half = 0.5 * rsq
#     ssq_half = 0.5 * ssq
#     rsoverssq = rsq / ssq
#     sigb0_inv_sigbq = solve(sigb0, sigbq)
#     -0.5 * (tr(WtW %*% sigbq) + tr(varphitvarphi %*% sigtq)) + sum(log(((pnorm(tmp))^y) * ((1 - pnorm(tmp))^(1-y)))) - 0.5 * J * (log(2 * pi) - (digamma(rsq_half) - log(ssq_half) + digamma(rtq_half) - log(stq_half))) + (0.25 * J * (J + 1) * ((sigpsiq * sqrt(2 / pi) * exp(-mupsiq^2/(2*sigpsiq2))) + (mupsiq * (1 - 2 * pnorm(-mupsiq / sigpsiq))))) - (0.5 * rsoverssq * rtoverstq * sum((diag(sigtq) + mutq^2) * (exp( (0.5 * sigpsiq2 * ((1:J)^2)) + (mupsiq * (1:J)) ) * (1 - pnorm( - mupsiq / sigpsiq - (sigpsiq * (1:J)))) + exp((0.5 * sigpsiq2 * ((1:J)^2)) - (mupsiq * (1:J))) * (1 - pnorm(mupsiq/sigpsiq - (sigpsiq * (1:J))))))) + (0.5 * J * (1 + log(2 * pi) + determinant(sigtq)$modulus[1])) + (rt0_half * log(st0_half)) - lgamma(rt0_half) + ((rt0_half + 1) * (digamma(rtq_half) - log(stq_half))) - st0_half * rtoverstq + rtq_half + log(stq_half) + lgamma(rtq_half) - ((1 + rtq_half) * digamma(rtq_half)) + rs0_half * log(ssq_half) - lgamma(rs0_half) + ((rs0_half + 1) * (digamma(rsq_half) - log(ssq_half))) - ss0_half * rsoverssq + rsq_half + log(ssq_half) + lgamma(rsq_half) - ((1 + rsq_half) * digamma(rsq_half)) + 0.5 * ( p + digamma(rsq_half) - log(ssq_half)  + determinant(sigb0_inv_sigbq)$modulus[1] - (rsoverssq * (tr(sigb0_inv_sigbq) + sum((mubq - mub0) * (solve(sigb0, mubq - mub0)))))) + log(0.5 * w0) - w0 * ((sigpsiq * sqrt(2 / pi) * exp(-mupsiq^2/(2 * sigpsiq2))) + (mupsiq * (1 - 2 * pnorm(-mupsiq/sigpsiq)))) + 0.5 * (log(2*pi*sigpsiq2) - 1)
#   }

#   VB = function(y, x, W, mub0, sigb0, w0, rs0, ss0, rt0, st0, mupsi0, J, tol = 1.0e-06) {
#     n = length(y)
#     if (!is.matrix(W)) W = as.matrix(W)
#     p = dim(W)[2]
#     dif = tol + 1
#     varphi = sqrt(2)*cos(outer(x,pi*(1:J)))
#     varphitvarphi = crossprod(varphi)
#     WtW = crossprod(W)
#     sigb0_inv = solve(sigb0)
#     sigb0_inv_mub0 = solve(sigb0, mub0)
#     # Initialize variational parameters
#     rsq = rs0 + J + p
#     ssq = ss0
#     rtq = rt0 + J
#     stq = st0
#     st0_half = st0/2
#     ss0_half = ss0/2
#     rt0_half = rt0/2
#     rs0_half = rs0/2
#     rtoverstq = rtq / stq
#     rsoverssq = rsq / ssq
#     mubq = mub0
#     mupsiq = mupsi0
#     sigpsiq2 = mupsi0^2 / 100
#     sigpsiq = sqrt(sigpsiq2)
#     mutq = rep(0, J)
#     muystar = rep(0, n)
#     sigbq = sigb0
#     sigtq = diag(1, J)
#     count = 0
#     lbold = LB(y, x, W, sigbq, varphi, sigtq, mubq, mutq, sigpsiq2, mupsiq, rsq, ssq, rtq, stq, rs0, ss0, rt0, st0, sigb0, mub0, w0)
#     lbnew = 0
#     lbrecord = c(lbold)
#     while (dif > tol) {
#       count = count + 1
#       a = 1
#       # Update psi
#       sigpsiq2_old = sigpsiq2
#       mupsiq_old = mupsiq
#       sigpsiq2_new = -0.5 / (der_S1_sigma2(sigpsiq2, sigpsiq, mupsiq, w0) + der_S2_sigma2(sigpsiq2, sigpsiq, sigtq, mutq, mupsiq, rsoverssq, rtoverstq, w0, J))
#       temp2 = sigpsiq2_new
#       mupsiq_old = mupsiq

#       # print(der_S2_mu(sigpsiq2, sigpsiq, sigtq, mutq, mupsiq, rsoverssq, rtoverstq, w0, J) == der_S1_mu(sigpsiq2, sigpsiq, mupsiq, w0))
#       temp = der_S1_mu(sigpsiq2, sigpsiq, mupsiq, w0) + der_S2_mu(sigpsiq2, sigpsiq, sigtq, mutq, mupsiq, rsoverssq, rtoverstq, w0, J)
#       while (sigpsiq2_new < 0) {
#         a = 2/3 * a
#         sigpsiq2_new = 1/(1/sigpsiq2_old + a * (1/sigpsiq2_new - 1/sigpsiq2_old))
#       }
#       sigpsiq2 = sigpsiq2_new
#       sigpsiq = sqrt(sigpsiq2)
#       mupsiq = mupsiq + a * sigpsiq2 * temp

#       # Update theta
#       sigtq = solve(varphitvarphi + (rsoverssq * rtoverstq * diag(MGF_foldedNormal(sigpsiq2, sigpsiq, mupsiq, (1:J)))))
#       mutq = sigtq %*% crossprod(varphi, muystar - W %*% mubq)

#       # Update tau^2
#       stq = st0 + (rtoverstq * (tr((sigtq + tcrossprod(mutq)) %*% diag(MGF_foldedNormal(sigpsiq2, sigpsiq, mupsiq, (1:J))))))
#       stq_half = stq / 2
#       rtoverstq = rtq / stq

#       # Update sigma^2
#       sigb0_inv_sigbq = solve(sigb0, sigbq)
#       ssq = ss0 + (rtoverstq * (tr((sigtq + tcrossprod(mutq)) %*% diag(MGF_foldedNormal(sigpsiq2, sigpsiq, mupsiq, (1:J)))))) + tr(sigb0_inv_sigbq) + sum((mubq - mub0) * solve(sigb0, (mubq-mub0)))
#       rsoverssq = rsq / ssq

#       # Update beta
#       sigbq = solve(WtW + rsoverssq * sigb0_inv)
#       mubq = sigbq %*% (rsoverssq * sigb0_inv_sigbq + crossprod(W, muystar - varphi %*% mutq))

#       # Update muystar
#       muystar_tmp = W %*% mubq + varphi %*% mutq
#       for (i in 1:n) {
#         if (y[i] == 1 & muystar_tmp[i] < -8) {
#           muystar[i] <- muystar_tmp[i] + dnormOverPnorm_Laurent(muystar_tmp[i])
#         } else if (y[i] == 1 & muystar_tmp[i] >= -8) {
#           muystar[i] <- muystar_tmp[i] + dnormOverPnorm(muystar_tmp[i])
#         } else {
#           muystar[i] <- muystar_tmp[i] + dnormOverPnormMinusOne(muystar_tmp[i])
#         }
#       }
#       # muystar = muystar_tmp + dnorm(muystar_tmp)/(((pnorm(muystar_tmp))^y)*((pnorm(muystar_tmp) - 1)^(1-y)))

#       lbnew = LB(y, x, W, sigbq, varphi, sigtq, mubq, mutq, sigpsiq2, mupsiq, rsq, ssq, rtq, stq, rs0, ss0, rt0, st0, sigb0, mub0, w0)
#       diff = lbnew - lbold
#       if (diff < 0) {
#         a = 1
#         sigpsiq2_new = temp2
#         while (sigpsiq2_new < 0) {
#           a = 2/3 * a
#           sigpsiq2_new = 1/(1/sigpsiq2 + a * (1/temp2 - 1/sigpsiq2))
#           # sigpsiq2_try = 1/(1/sigpsiq2_old + step * (1/sigpsiq2 - 1/sigpsiq2_old))
#           # sigpsiq_try = sqrt(sigpsiq2_try)
#           # mupsiq_try = sigpsiq2_try * (mupsiq_old/sigpsiq2_old + step*(mupsiq/sigpsiq2 - mupsiq_old/sigpsiq2_old))
#           # lbnew = LB(y, x, W, sigbq, varphi, sigtq, mubq, mutq, sigpsiq2_try, mupsiq_try, rsq, ssq, rtq, stq, rs0, ss0, rt0, st0, sigb0, mub0, w0)
#           # dif_try = lbnew - lbtmp
#           # cat('sigpsiq2_try: ', sigpsiq2_try, ', mupsiq_try: ', mupsiq_try, '\n')
#         }
#         sigpsiq2 = sigpsiq2_new
#         sigpsiq = sqrt(sigpsiq2)
#         mupsiq = mupsiq_old + a * sigpsiq2 * temp

#         # Update theta
#         sigtq = solve(varphitvarphi + (rsoverssq * rtoverstq * diag(MGF_foldedNormal(sigpsiq2, sigpsiq, mupsiq, (1:J)))))
#         mutq = sigtq %*% crossprod(varphi, muystar - W %*% mubq)

#         # Update tau^2
#         stq = st0 + (rtoverstq * (tr((sigtq + tcrossprod(mutq)) %*% diag(MGF_foldedNormal(sigpsiq2, sigpsiq, mupsiq, (1:J))))))
#         stq_half = stq / 2
#         rtoverstq = rtq / stq

#         # Update sigma^2
#         sigb0_inv_sigbq = solve(sigb0, sigbq)
#         ssq = ss0 + (rtoverstq * (tr((sigtq + tcrossprod(mutq)) %*% diag(MGF_foldedNormal(sigpsiq2, sigpsiq, mupsiq, (1:J)))))) + tr(sigb0_inv_sigbq) + sum((mubq - mub0) * solve(sigb0, (mubq-mub0)))
#         rsoverssq = rsq / ssq

#         # Update beta
#         sigbq = solve(rsoverssq * (WtW + sigb0_inv))
#         mubq = sigbq %*% (rsoverssq * sigb0_inv_sigbq + crossprod(W, muystar - varphi %*% mutq))

#         # Update muystar
#         muystar_tmp = W %*% mubq + varphi %*% mutq
#         for (i in 1:n) {
#         if (y[i] == 1 & muystar_tmp[i] < -8) {
#           muystar[i] <- muystar_tmp[i] + dnormOverPnorm_Laurent(muystar_tmp[i])
#         } else if (y[i] == 1 & muystar_tmp[i] >= -8) {
#           muystar[i] <- muystar_tmp[i] + dnormOverPnorm(muystar_tmp[i])
#         } else {
#           muystar[i] <- muystar_tmp[i] + dnormOverPnormMinusOne(muystar_tmp[i])
#         }
#       }
#         # muystar = muystar_tmp + dnorm(muystar_tmp)/(((pnorm(muystar_tmp))^y)*((pnorm(muystar_tmp) - 1)^(1-y)))

#         lbnew = LB(y, x, W, sigbq, varphi, sigtq, mubq, mutq, sigpsiq2, mupsiq, rsq, ssq, rtq, stq, rs0, ss0, rt0, st0, sigb0, mub0, w0)
#         diff = lbnew - lbold
#       }

#       dif = (lbnew - lbold)/abs(lbnew)
#       lbold = lbnew
#       lbrecord = c(lbrecord, lbnew)
#       cat("count: ", count, ", lbnew: ", lbnew, ", dif: ", dif, "\n")
#     }
#     list(mutq = mutq, sigtq = sigtq, mubq = mubq, sigbq = sigbq, rtq = rtq, stq = stq, rsq = rsq, ssq = ssq, muystar = muystar, lbrecord = lbrecord)
#   }


#   ##########################################################

#   ################    Start simulation    ##################

#   ##########################################################

#   N_ = 1000
#   x = runif(N_, 0, 1)
#   W = rep(1, times = N_)
#   mub0 = 0
#   y = rbinom(length(x),1,pnorm(FUN(x) + W))
#   sigb0 = matrix(1,nrow=1,ncol=1)
#   fit = VB(y, x, W, mub0, sigb0, w0 = 1, rs0 = 0.01, ss0 = 0.01, rt0 = 0.01, st0 = 0.01, mupsi0 = 1, J = J, tol = 1.0e-06)
#   categorize = function(x) {
#     y = 0
#     if (x < 0) {
#       y = 0
#     } else {
#       y = 1
#     }
#     y
#   }

#   varphi = sqrt(2)*cos(outer(x,pi*(1:J)))
#   unknown_f = varphi %*% fit$mutq + W %*% fit$mubq
#   res = sapply(fit$muystar, categorize)
#   if (draw == TRUE) {
#     # plot(y-mean(y) ~ xobs, xlab = 'index', ylab = 'observed/fitted', main = 'Simulation result', type = 'p')
#     ord = order(x)
#     plot(x[ord], unknown_f[ord], col = 'purple', xlab = 'index', ylab = 'observed/fitted', main = 'Simulation result', ylim = c(-10, 10), type = 'l')
#     # lines(xobs[ord], fixed[ord], col = 'purple')
#     curve(FUN, from = 0, to = 1, col = 'red', lty = 2, add = TRUE)
#     legend("topright", legend = c('fitted values', 'true function'), col = c('purple', 'red'), lty = c(1, 2), bg = 'gray95')
#     return(list(fit = fit, res = res, y = y, x = x))
#   } else {
#     return(list(fit = fit, res = res, y = y, x = x))
#   }
# }

# # sim_cosine_probit = function(FUN, J) {
# #   der_Qj_mu = function(sigma2, sigma, mu, j) {
# #     (exp((0.5 * sigma2 * j^2) + (mu * j)) * dnorm(-mu/sigma - sigma * j) / sigma) + (j * exp((0.5 * sigma2 * j^2) + (mu * j)) * (1 - pnorm((-mu/sigma) - (sigma * j)))) - (exp(0.5 * sigma2 * j^2 - mu * j) * dnorm(mu/sigma - sigma * j) / sigma) - (j * exp((0.5 * sigma2 * j^2) - (mu * j)) * (1 - pnorm(mu/sigma - sigma * j)))
# #   }
# #   der_Qj_sigma2 = function(sigma2, sigma, mu, j) {
# #     (-exp(0.5 * sigma2 * j^2 + mu * j) * dnorm(-mu/sigma - sigma * j) * ((mu/(2 * sigma^3)) - (j / (2 * sigma2)))) + (0.5 * j^2 * exp(0.5 * sigma2 * j^2 + mu * j) * (1 - pnorm(-mu/sigma - (sigma * j)))) + (exp((0.5 * sigma2 * j^2) - (mu * j)) * dnorm(mu/sigma - (sigma * j)) * ((mu/(2 * sigma^3)) + (j/(2 * sigma2)))) + (0.5 * j^2 * exp((0.5 * sigma2 * j^2) - (mu * j)) * (1 - pnorm(mu/sigma - sigma * j)))
# #   }

# #   der_S1_mu = function(sigma2, sigma, mu, w0) {
# #     -w0 * ( (-mu/sigma * sqrt(2/pi) * exp(-mu^2 / (2 * sigma2))) + (1 - 2 * pnorm(-mu/sigma)) + (2*mu/sigma * dnorm(-mu/sigma)) )
# #   }

# #   der_S1_sigma2 = function(sigma2, sigma, mu, w0) {
# #     -w0 * (((1/(2 * sigma)) + (0.5 * mu^2 / sigma^3) * sqrt(2/pi) * exp(-mu^2 / (2 * sigma2))) - (mu^2 / sigma^3 * dnorm(-mu/sigma)))
# #   }

# #   der_S2_mu = function(sigma2, sigma, sigtq, mutq, mu, rsoverssq, rtoverstq, w0, J) {
# #     ((-0.25 * J * (J + 1) / w0) * der_S1_mu(sigma2, sigma, mu, w0)) - (0.5 * rsoverssq * rtoverstq * sum(((diag(sigtq) + mutq^2) * der_Qj_mu(sigma2, sigma, mu, (1:J)))))
# #   }

# #   der_S2_sigma2 = function(sigma2, sigma, sigtq, mutq, mu, rsoverssq, rtoverstq, w0, J) {
# #     ((-0.25 * J * (J + 1) / w0) * der_S1_sigma2(sigma2, sigma, mu, w0)) - (0.5 * rsoverssq * rtoverstq * sum(((diag(sigtq) + mutq^2) * der_Qj_sigma2(sigma2, sigma, mu, (1:J)))))
# #   }
# #   Qj = function(mu, sigma2, j) {
# #     term1 = 0.5 * sigma2 * j^2
# #     term2 = mu * j
# #     term3 = mu/sqrt(sigma2)
# #     term4 = sqrt(sigma2) * j
# #     (exp(term1 + term2) * (1 - pnorm(-(term3 + term4)))) + (exp(term1 - term2) * (1 - pnorm(term3 - term4)))
# #   }
# #   S1 = function(mu, sigma2, w0) -w0 * ((sqrt(sigma2 * 2 / pi) * exp(-mu^2 / (2 * sigma2))) + (mu * (1 - 2 * pnorm(-mu/sqrt(sigma2)))))
# #   S2 = function(mupsiq, sigpsiq2, w0, J, rtq, stq, rsq, ssq, sigtq, mutq) ((-0.5 * rtq/stq * rsq/ssq) * sum((diag(sigtq) + mutq^2) * Qj(mupsiq, sigpsiq2, 1:J))) - (J*(J+1)/(4*w0) * S1(mupsiq, sigpsiq2, w0))

# #   LB = function(y, W, varphi, mubq, mutq, sigb0, sigbq, rsq, ssq, rs0, ss0, rtq, stq, rt0, st0, mub0, mupsiq, sigpsiq2, w0, J) {
# #     p = dim(W)[2]
# #     st0_half = st0/2
# #     ss0_half = ss0/2
# #     rt0_half = rt0/2
# #     rs0_half = rs0/2
# #     stq_half = stq/2
# #     ssq_half = ssq/2
# #     rtq_half = rtq/2
# #     rsq_half = rsq/2
# #     tau_ratio0 = rt0/st0
# #     sig_ratio0 = rs0/ss0
# #     tau_ratioq = rtq/stq
# #     sig_ratioq = rsq/ssq
# #     WtW = crossprod(W)
# #     PtP = crossprod(varphi)
# #     t1 = W %*% mubq + varphi %*% mutq
# #     t2 = solve(sigb0, sigbq)
# #     (-0.5 * sum(WtW * sigbq) + sum(PtP * sigtq)) + sum(log((pnorm(t1)^y) * ((1-pnorm(t1))^(1-y))))
# #     + 0.5 * (J * (digamma(rsq_half) - log(ssq_half) + digamma(rtq_half) - log(stq_half) + 1) + determinant(sigtq, logarithm = TRUE)$modulus[1] + determinant(t2, logarithm = TRUE)$modulus[1])
# #     + rt0_half * log(st0_half) - lgamma(rt0_half) + (rt0_half + 1) * (digamma(rtq_half) - log(stq_half)) - st0_half * tau_ratioq + rtq_half + log(stq_half) + lgamma(rtq_half)
# #     - (1 + rtq_half) * digamma(rtq_half)
# #     +rs0_half * log(ss0_half) - lgamma(rs0_half) + (rs0_half + 1) * (digamma(rsq_half)-log(ssq_half)) - ss0_half * sig_ratioq + rsq_half +log(ssq_half) + lgamma(rsq_half)
# #     - (1 + rsq_half) * digamma(rsq_half)
# #     + 0.5 * (p + digamma(rsq_half) - log(ssq_half) - sig_ratioq * (sum(diag(t2)) + sum((mubq - mub0) * solve(sigb0, mubq - mub0))) + log(2*pi*sigpsiq2) + 1) + log(w0/2)
# #     + S1(mupsiq, sigpsiq2, w0) + S2(mupsiq, sigpsiq2, w0, J, rtq, stq, rsq, ssq, sigtq, mutq)
# #   }

# #   VB = function(y, W, x, st0, st0, rs0, ss0, mub0, sigb0, w0, J, tol = 1.0e-6) {
# #     WtW = crossprod(W)
# #     PtP = crossprod(varphi)
# #     p = dim(W)[2]
# #     # Initialize variational parameters
# #     rtq = rt0 + J
# #     stq = 0.01
# #     tau_ratioq = rtq/stq
# #     rsq = rs0 + J + p
# #     ssq = 0.01
# #     sig_ratioq = rsq/ssq
# #     mutq = rep(0, J)
# #     sigtq = diag(1, J)
# #     sigbq = diag(1, p)
# #     mubq = rep(0, p)
# #     sigpsiq2 = 0.01
# #     mupsiq = 0.1
# #     muystar = rep(0, length(y))

# #     sigb0_inv = solve(sigb0)

# #     dif = tol + 1
# #     lbold = 10e-7
# #     lbnew = 0
# #     lbrecord = c()
# #     count = 0
# #     while (dif > tol) {
# #       count = count + 1
# #       a = 1
# #       # Update psi
# #       sigpsiq2_old = sigpsiq2
# #       sigpsiq2_new = -0.5 / ((der_S1_sigma2(sigpsiq2, sigpsiq, mupsiq, w0) + der_S2_sigma2(sigpsiq2, sigpsiq, sigtq, mutq, mupsiq, sig_ratioq, tau_ratioq, w0, J)))
# #       temp2 = sigpsiq2_new
# #       mupsiq_old = mupsiq
# #       temp = der_S1_mu(sigpsiq2, sigpsiq, mupsiq, w0) + der_S2_mu(sigpsiq2, sigpsiq, sigtq, mutq, mupsiq, sig_ratioq, tau_ratioq, w0, J)
# #       while (sigpsi2_new < 0) {
# #         a = 2/3 * a
# #         sigpsiq2_new = 1/(1/sigpsiq2_old + a * (1/sigpsiq2_new - 1/sigpsiq2_old))
# #       }
# #       sigpsiq2 = sigpsiq2_new
# #       sigpsiq = sqrt(sigpsiq2)
# #       mupsiq = mupsiq + a * sigpsiq2 * temp

# #       # Update theta
# #       sigtq = solve(PtP + sig_ratioq * tau_ratioq * diag(Qj(mupsiq, sigpsiq2, 1:J)))
# #       mutq = sigtq %*% crossprod(varphi, muystar - W %*% mubq)
      
# #       # Update tau
# #       stq = st0 + sig_ratioq * sum((diag(sigtq) + mutq^2) * Qj(mupsiq, sigpsiq2, 1:J))
# #       tau_ratioq = rtq/stq

# #       # Update sigma
# #       ssq = ss0 + tau_ratioq * sum((diag(sigtq) + mutq^2) * Qj(mupsiq, sigpsiq2, 1:J)) + sum(sigb0_inv * sigbq) + sum((mubq - mub0) * (sigb0_inv %*% (mubq - mub0)))
# #       sig_ratioq = rsq/ssq
      
# #       # Update beta
# #       sigbq = solve(sig_ratioq * (Wtw + sigb0_inv))
# #       mubq = sig_ratioq * sigbq %*% (sigb0_inv %*% mub0 + crossprod(W, muystar - varphi %*% mutq))

# #       # Update muystar

# #     }


# #   }
# # }


# require(Rcpp)
# require(RcppArmadillo)
# sourceCpp('auxiliary.cpp')

# dnormOverPnorm <- function(x) sqrt(2 / pi) * exp(-x^2 / 2) / (erf(x / sqrt(2)) + 1)
# dnormOverPnorm_Laurent <- function(x) -x - 1 / x + 2 / x^3 - 10 / x^5
# dnormOverPnormMinusOne <- function(x) exp(-x^2 / 2) / (sqrt(pi / 2) * (erf(x / sqrt(2)) + 1) - 1)

Qj <- function(muPsi, sigmaPsi2, sigmaPsi, j) {
  t1 <- 0.5 * sigmaPsi2 * j^2
  t2 <- muPsi * j
  t3 <- muPsi / sigmaPsi
  t4 <- sigmaPsi * j
  exp(t1 + t2) * (1 - pnorm(-t3 - t4)) + exp(t1 - t2) * (1 - pnorm(t3 - t4))
}

S1 <- function(muPsi, sigmaPsi2, sigmaPsi, w0) {
  -w0 * (sigmaPsi * sqrt(2/pi) * exp(-0.5 * muPsi^2 / sigmaPsi2) + muPsi * (1 - 2 * pnorm(-muPsi / sigmaPsi)))
}

S2 <- function(muPsi, sigmaPsi2, sigmaPsi, w0, rqs, sqs, rqt, sqt, SigmaThetaq, muThetaq, J) {
  -0.5 * rqs/sqs * rqt/sqt * sum((diag(SigmaThetaq) + muThetaq^2) * Qj(muPsi, sigmaPsi2, sigmaPsi, 1:J)) - 0.25 * J * (J + 1) / w0 * S1(muPsi, sigmaPsi2, sigmaPsi, w0)
}

DQj_muPsi <- function(muPsi, sigmaPsi2, sigmaPsi, j) {
  t1 <- 0.5 * sigmaPsi2 * j^2
  t2 <- muPsi * j
  t3 <- muPsi / sigmaPsi
  t4 <- sigmaPsi * j
  t5 <- 1 / sigmaPsi
  exp(t1 + t2) * dnorm(-t3 - t4) * t5 + j * exp(t1 + t2) * (1 - pnorm(-t3 - t4)) - exp(t1 - t2) * dnorm(t3 - t4) * t5 - j * exp(t1 - t2) * (1 - pnorm(t3 - t4))
}

DQj_sigmaPsi2 <- function(muPsi, sigmaPsi2, sigmaPsi, j) {
  t1 <- 0.5 * sigmaPsi2 * j^2
  t2 <- muPsi * j
  t3 <- muPsi / sigmaPsi
  t4 <- sigmaPsi * j
  t5 <- 0.5 * muPsi / sigmaPsi^3
  t6 <- 0.5 * j / sigmaPsi
  t7 <- 0.5 * j^2
  -exp(t1 + t2) * dnorm(-t3 - t4) * (t5 - t6) + t7 * exp(t1 + t2) * (1 - pnorm(-t3 - t4)) + exp(t1 - t2) * dnorm(t3 - t4) * (t5 + t6) + t7 * exp(t1 - t2) * (1 - pnorm(t3 - t4))
}

DS1_muPsi <- function(muPsi, sigmaPsi2, sigmaPsi, w0) {
  -w0 * ((0.5 / sigmaPsi * (1 + muPsi^2 / sigmaPsi2)) * sqrt(2 / pi) * exp(-0.5 * muPsi^2 / sigmaPsi2) - muPsi^2 / sigmaPsi^3 * dnorm(-muPsi / sigmaPsi))
}

DS1_sigmaPsi2 <- function(muPsi, sigmaPsi2, sigmaPsi, w0) {
  -w0 * (-muPsi / sigmaPsi * sqrt(2 / pi) * exp(-0.5 * muPsi^2 / sigmaPsi2) + (1 - 2 * pnorm(-muPsi / sigmaPsi)) + 2 * muPsi / sigmaPsi * dnorm(-muPsi / sigmaPsi))
}

DS2_muPsi <- function(muPsi, sigmaPsi2, sigmaPsi, w0, rqs, sqs, rqt, sqt, SigmaThetaq, muThetaq, J) {
  -0.25 * J * (J + 1) / w0 * DS1_sigmaPsi2(muPsi, sigmaPsi2, sigmaPsi, w0) - 0.5 * rqs / sqs * rqt / sqt * sum((diag(SigmaThetaq) + muThetaq^2) * DQj_muPsi(muPsi, sigmaPsi2, sigmaPsi, 1:J))
}

DS2_sigmaPsi2 <- function(muPsi, sigmaPsi2, sigmaPsi, w0, rqs, sqs, rqt, sqt, SigmaThetaq, muThetaq, J) {
  -0.25 * J * (J + 1) / w0 * DS1_sigmaPsi2(muPsi, sigmaPsi2, sigmaPsi, w0) - 0.5 * rqs / sqs * rqt / sqs * sum((diag(SigmaThetaq) + muThetaq^2) * DQj_sigmaPsi2(muPsi, sigmaPsi2, sigmaPsi, 1:J))
}

LB_experiment <- function(y, W, varphi, WtW, varphitvarphi, muPsi, sigmaPsi2, sigmaPsi, w0, J, SigmaThetaq, muThetaq, SigmaBetaq, muBetaq, SigmaBeta0_inv, muBeta0, rqt, sqt, rqt_half, sqt_half, rqs, sqs, rqs_half, sqs_half, r0t_half, s0t_half, r0s_half, s0s_half) {
  p <- dim(W)[2]
  t1 <- SigmaBeta0_inv %*% SigmaBetaq
  -0.5 * (sum(WtW * SigmaBetaq) + sum(varphitvarphi * SigmaThetaq)) + sum(log(pnorm(W %*% muBetaq + varphi %*% muThetaq)^y * (1 - pnorm(W %*% muBetaq + varphi %*% muThetaq))^(1 - y))) - 0.5 * J * (log(2 * pi) - (digamma(rqs_half) - log(sqs_half)) - (digamma(rqt_half) - log(sqt_half))) + S2(muPsi, sigmaPsi2, sigmaPsi, w0, rqs, sqs, rqt, sqt, SigmaThetaq, muThetaq, J) + 0.5 * J * (1 + log(2 * pi)) + 0.5 * determinant(SigmaThetaq)$modulus[1] + r0t_half * log(r0t_half) - lgamma(r0t_half) + (r0t_half - rqt_half) * digamma(rqt_half) - (r0t_half + 1) * log(sqt_half) - r0t_half * rqt / sqt + rqt_half + log(sqt_half) + lgamma(rqt_half) + r0s_half * log(s0s_half) - lgamma(r0s_half) + (r0s_half - rqs_half) * digamma(rqs_half) - (r0s_half + 1) * log(sqs_half) - sqs_half * rqs / sqs + rqs_half + log(sqs_half) + lgamma(rqs_half) + 0.5 * (p + 1) * (1 + digamma(rqs_half) - log(sqs_half)) + 0.5 * determinant(t1)$modulus[1] - 0.5 * rqs / sqs * (sum(diag(t1)) + sum((muBetaq - muBeta0) * (SigmaBeta0_inv %*% (muBetaq - muBeta0)))) + log(w0 / 2) + S1(muPsi, sigmaPsi2, sigmaPsi, w0) + 0.5 * (log(2 * pi * sigmaPsi2) - 1)
}

VB_experiment <- function(x, y, W, muBeta0, SigmaBeta0, w0, r0s, s0s, r0t, s0t, J, tol = 1.0e-06) {
  n <- length(y)
  if (!is.matrix(W)) W <- as.matrix(W)
  p <- dim(W)[2]
  dif <- tol + 1
  varphi <- sqrt(2) * cos(outer(x, pi * (1:J)))
  varphitvarphi <- crossprod(varphi)
  WtW <- crossprod(W)
  SigmaBeta0_inv <- solve(SigmaBeta0)
  SigmaBeta0_inv_muBeta0 <- solve(SigmaBeta0, muBeta0)
  rqs <- r0s + J + p
  rqs_half <- 0.5 * rqs
  sqs <- s0s + 100
  sqs_half <- 0.5 * sqs
  rqt <- r0t + J
  rqt_half <- 0.5 * rqt
  sqt <- s0t + 100
  sqt_half <- 0.5 * sqt
  s0t_half <- 0.5 * s0t
  s0s_half <- 0.5 * s0s
  r0t_half <- 0.5 * r0t
  r0s_half <- 0.5 * r0s
  muBetaq <- muBeta0
  muPsi <- 0.1
  sigmaPsi2 <- muPsi^2 / 100
  sigmaPsi <- sqrt(sigmaPsi2)
  muThetaq <- rep(0, J)
  muYStar <- rep(0, n)
  muYStarq <- rep(1, n)
  SigmaBetaq <- SigmaBeta0
  SigmaThetaq <- diag(1, J)
  count <- 0
  lbold <- LB_experiment(y, W, varphi, WtW, varphitvarphi, muPsi, sigmaPsi2, sigmaPsi, w0, J, SigmaThetaq, muThetaq, SigmaBetaq, muBetaq, SigmaBeta0_inv, muBeta0, rqt, sqt, rqt_half, sqt_half, rqs, sqs, rqs_half, sqs_half, r0t_half, s0t_half, r0s_half, s0s_half)
  cat('first lb: ', lbold, '\n')
  lbnew <- 0
  lbrecord <- c(lbold)
  while (dif > tol | dif < 0) {
    count <- count + 1

    # Update theta
    SigmaThetaq <- solve(varphitvarphi + rqs / sqs * rqt / sqt * diag(Qj(muPsi, sigmaPsi2, sigmaPsi, 1:J)))
    muThetaq <- SigmaThetaq %*% crossprod(varphi, muYStarq - W %*% muBetaq)

    # Update tau
    sqt <- s0t + rqs / sqs * sum((diag(SigmaThetaq) + muThetaq^2) * Qj(muPsi, sigmaPsi2, sigmaPsi, 1:J))
    sqt_half <- 0.5 * sqt

    # Update sigma
    cat('sum((diag(SigmaThetaq) + muThetaq^2) * Qj(muPsi, sigmaPsi2, sigmaPsi, 1:J)):', sum((diag(SigmaThetaq) + muThetaq^2) * Qj(muPsi, sigmaPsi2, sigmaPsi, 1:J)), '\n')
    cat('sum(SigmaBeta0_inv * SigmaBetaq): ', sum(SigmaBeta0_inv * SigmaBetaq), '\n')
    sqs <- s0s + rqt / sqt * sum((diag(SigmaThetaq) + muThetaq^2) * Qj(muPsi, sigmaPsi2, sigmaPsi, 1:J)) + sum(SigmaBeta0_inv * SigmaBetaq) + sum((muBetaq - muBeta0) * (SigmaBeta0_inv %*% (muBetaq - muBeta0)))

    # Update beta
    cat('rqs/sqs: ', rqs/sqs, '\n')
    cat('SigmaBeta0_inv: ', SigmaBeta0_inv, '\n')
    cat('WtW: ', WtW, '\n')
    cat('rqs/sqs * SigmaBeta0_inv + WtW: ', rqs / sqs * SigmaBeta0_inv + WtW, '\n')
    SigmaBetaq <- solve(as.matrix(rqs / sqs * SigmaBeta0_inv + WtW))
    muBetaq <- SigmaBetaq %*% (rqs / sqs * SigmaBeta0_inv_muBeta0 + crossprod(muYStar - varphi %*% muThetaq))

    # Update muYStar
    muYStar <- W %*% muBetaq + varphi %*% muThetaq
    # for (i in 1:n) {
    #   if (y[i] == 1 & muYStar[i] < -8) {
    #     muYStarq[i] <- muYStar[i] + dnormOverPnorm_Laurent(muYStar[i])
    #   } else if (y[i] == 1 & muYStar[i] >= -8) {
    #     muYStarq[i] <- muYStar[i] + dnormOverPnorm(muYStar[i])
    #   } else {
    #     muYStarq[i] <- muYStar[i] + dnormOverPnormMinusOne(muYStar[i])
    #   }
    # }
    muYStarq <- muYStar + dnorm(muYStar)/((pnorm(muYStar)^y)*(pnorm(muYStar)-1)^(1-y))
    a <- 1
    # Update psi
    sigmaPsi2_old <- sigmaPsi2
    muPsi_old <- muPsi
    temp <- DS1_muPsi(muPsi, sigmaPsi2, sigmaPsi, w0) + DS2_muPsi(muPsi, sigmaPsi2, sigmaPsi, w0, rqs, sqs, rqt, sqt, SigmaThetaq, muThetaq, J)
    # sigmaPsi2 <- -0.5 / (DS1_sigmaPsi2(muPsi, sigmaPsi2, sigmaPsi, w0) + DS2_sigmaPsi2(muPsi, sigmaPsi2, sigmaPsi, w0, rqs, sqs, rqt, sqt, SigmaThetaq, muThetaq, J))
    
    sigmaPsi2_new <- -0.5 / (DS1_sigmaPsi2(muPsi, sigmaPsi2, sigmaPsi, w0) + DS2_sigmaPsi2(muPsi, sigmaPsi2, sigmaPsi, w0, rqs, sqs, rqt, sqt, SigmaThetaq, muThetaq, J))
    # temp2 <- sigmaPsi2_new

    while (sigmaPsi2_new < 0) {
      a <- 2/3 * a
      sigmaPsi2_new <- 1 / (1 / sigmaPsi2_old + a * (1 / sigmaPsi2_new - 1 / sigmaPsi2_old))
    }
    sigmaPsi2 <- sigmaPsi2_new
    sigmaPsi <- sqrt(sigmaPsi2)
    muPsi <- muPsi + a * sigmaPsi2 * temp
    cat('sigmaPsi2: ', sigmaPsi2, '\n')
    cat('muPsi: ', muPsi, '\n')
    lbnew <- LB_experiment(y, W, varphi, WtW, varphitvarphi, muPsi, sigmaPsi2, sigmaPsi, w0, J, SigmaThetaq, muThetaq, SigmaBetaq, muBetaq, SigmaBeta0_inv, muBeta0, rqt, sqt, rqt_half, sqt_half, rqs, sqs, rqs_half, sqs_half, r0t_half, s0t_half, r0s_half, s0s_half)
    lbfull <- lbnew
    dif <- lbnew - lbold
    cat('lbfull: ', lbfull, '\n')

    if (dif < 0) {
      a <- 1
      dif_try <- dif
      while (dif_try < 0) {
        a <- 0.5 * a
        sigmaPsi2_try <- 1 / (1 / sigmaPsi2_old + a * (1 / sigmaPsi2 - 1 / sigmaPsi2_old))
        sigmaPsi_try <- sqrt(sigmaPsi2_try)
        muPsi_try <- sigmaPsi2_try * (muPsi_old / sigmaPsi2_old + a * (muPsi / sigmaPsi2 - muPsi_old / sigmaPsi2_old))

        # Update theta
        SigmaThetaq <- solve(varphitvarphi + rqs / sqs * rqt / sqt * diag(Qj(muPsi, sigmaPsi2, sigmaPsi, 1:J)))
        muThetaq <- SigmaThetaq %*% crossprod(varphi, muYStarq - W %*% muBetaq)

        # Update tau
        sqt <- s0t + rqs / sqs * sum((diag(SigmaThetaq) + muThetaq^2) * Qj(muPsi, sigmaPsi2, sigmaPsi, 1:J))
        sqt_half <- 0.5 * sqt

        # Update sigma
        cat('sum((diag(SigmaThetaq) + muThetaq^2) * Qj(muPsi, sigmaPsi2, sigmaPsi, 1:J)):', sum((diag(SigmaThetaq) + muThetaq^2) * Qj(muPsi, sigmaPsi2, sigmaPsi, 1:J)), '\n')
        cat('sum(SigmaBeta0_inv * SigmaBetaq): ', sum(SigmaBeta0_inv * SigmaBetaq), '\n')
        sqs <- s0s + rqt / sqt * sum((diag(SigmaThetaq) + muThetaq^2) * Qj(muPsi, sigmaPsi2, sigmaPsi, 1:J)) + sum(SigmaBeta0_inv * SigmaBetaq) + sum((muBetaq - muBeta0) * (SigmaBeta0_inv %*% (muBetaq - muBeta0)))

        # Update beta
        cat('rqs/sqs: ', rqs/sqs, '\n')
        cat('SigmaBeta0_inv: ', SigmaBeta0_inv, '\n')
        cat('WtW: ', WtW, '\n')
        cat('rqs/sqs * SigmaBeta0_inv + WtW: ', rqs / sqs * SigmaBeta0_inv + WtW, '\n')
        SigmaBetaq <- solve(as.matrix(rqs / sqs * SigmaBeta0_inv + WtW))
        muBetaq <- SigmaBetaq %*% (rqs / sqs * SigmaBeta0_inv_muBeta0 + crossprod(muYStar - varphi %*% muThetaq))

        # Update muYStar
        muYStar <- W %*% muBetaq + varphi %*% muThetaq
        # for (i in 1:n) {
        #   if (y[i] == 1 & muYStar[i] < -8) {
        #     muYStarq[i] <- muYStar[i] + dnormOverPnorm_Laurent(muYStar[i])
        #   } else if (y[i] == 1 & muYStar[i] >= -8) {
        #     muYStarq[i] <- muYStar[i] + dnormOverPnorm(muYStar[i])
        #   } else {
        #     muYStarq[i] <- muYStar[i] + dnormOverPnormMinusOne(muYStar[i])
        #   }
        # }
        muYStarq <- muYStar + dnorm(muYStar)/((pnorm(muYStar)^y)*(pnorm(muYStar)-1)^(1-y))


        lbnew <- LB_experiment(y, W, varphi, WtW, varphitvarphi, muPsi, sigmaPsi2, sigmaPsi, w0, J, SigmaThetaq, muThetaq, SigmaBetaq, muBetaq, SigmaBeta0_inv, muBeta0, rqt, sqt, rqt_half, sqt_half, rqs, sqs, rqs_half, sqs_half, r0t_half, s0t_half, r0s_half, s0s_half)
        dif_try <- lbnew - lbfull
      }

      sigmaPsi2 <- sigmaPsi2_try
      sigmaPsi <- sigmaPsi_try
      muPsi <- muPsi_try

    #   # Update psi
    #   sigmaPsi2_old <- sigmaPsi2
    #   muPsi_old <- muPsi
    #   sigmaPsi2_new <- -0.5 / (DS1_sigmaPsi2(muPsi, sigmaPsi2, sigmaPsi, w0) + DS2_sigmaPsi2(muPsi, sigmaPsi2, sigmaPsi, w0, rqs, sqs, rqt, sqt, SigmaThetaq, muThetaq, J))
    #   temp2 <- sigmaPsi2_new

    #   temp <- DS1_muPsi(muPsi, sigmaPsi2, sigmaPsi, w0) + DS2_muPsi(muPsi, sigmaPsi2, sigmaPsi, w0, rqs, sqs, rqt, sqt, SigmaThetaq, muThetaq, J)
    #   while (sigmaPsi2_new < 0) {
    #     a <- 2/3 * a
    #     sigmaPsi2_new <- 1 / (1 / sigmaPsi2_old + a * (1 / sigmaPsi2_new - 1 / sigmaPsi2_old))
    #   }
    #   sigmaPsi2 <- sigmaPsi2_new
    #   sigmaPsi <- sqrt(sigmaPsi2)
    #   muPsi <- muPsi + a * sigmaPsi2 * temp

    #   # Update theta
    #   SigmaThetaq <- solve(varphitvarphi + rqs / sqs * rqt / sqt * diag(Qj(muPsi, sigmaPsi2, sigmaPsi, 1:J)))
    #   muThetaq <- SigmaThetaq %*% crossprod(varphi, muYStarq - W %*% muBetaq)

    #   # Update tau
    #   sqt <- s0t + rqs / sqs * sum((diag(SigmaThetaq) + muThetaq^2) * Qj(muPsi, sigmaPsi2, sigmaPsi, 1:J))
    #   sqt_half <- 0.5 * sqt

    #   # Update sigma
    #   sqs <- s0s + rqt / sqt * sum((diag(SigmaThetaq) + muThetaq^2) * Qj(muPsi, sigmaPsi2, sigmaPsi, 1:J)) + sum(SigmaBeta0_inv * SigmaBetaq) + sum((muBetaq - muBeta0) * (SigmaBeta0_inv %*% (muBetaq - muBeta0)))

    #   # Update beta
    #   SigmaBetaq <- solve(rqs / sqs * SigmaBeta0_inv + WtW)
    #   muBetaq <- SigmaBetaq %*% (rqs / sqs * SigmaBeta0_inv_muBeta0 + crossprod(muYStar - varphi %*% muThetaq))

    #   # Update muYStar
    #   muYStar <- W %*% muBetaq + varphi %*% muThetaq
    #   for (i in 1:n) {
    #     if (y[i] > 0 & muYStar[i] < -8) {
    #       muYStarq[i] <- muYStar[i] + dnormOverPnorm_Laurent(muYStar[i])
    #     } else if (y[i] > 0 & muYStar[i] >= -8) {
    #       muYStarq[i] <- muYStar[i] + dnormOverPnorm(muYStar[i])
    #     } else {
    #       muYStarq[i] <- muYStar[i] + dnormOverPnormMinusOne(muYStar[i])
    #     }
    #   }
    #   lbnew <- LB(y, W, varphi, WtW, varphitvarphi, muPsi, sigmaPsi2, sigmaPsi, w0, J, SigmaThetaq, muThetaq, SigmaBetaq, muBetaq, SigmaBeta0_inv, muBeta0, rqt, sqt, rqt_half, sqt_half, rqs, sqs, rqs_half, sqs_half, r0t_half, s0t_half, r0s_half, s0s_half)
    #   diff <- lbnew - lbold
    # }

    # dif <- (lbnew - lbold) / abs(lbnew)
    # lbold <- lbnew
    # lbrecord <- c(lbrecord, lbnew)
    }
    dif <- (lbnew - lbold)/abs(lbnew)
    lbold <- lbnew
    lbrecord <- c(lbrecord, lbnew)
    cat('count: ', count, ', lbnew: ', lbnew, ', dif: ', dif, '\n')
  }
  list(muThetaq = muThetaq, SigmaThetaq = SigmaThetaq, muBetaq = muBetaq, SigmaBetaq = SigmaBetaq, rqt = rqt, sqt = sqt, rqs = rqs, sqs = sqs, muYStarq = muYStarq, lbrecord = lbrecord)
}

simExperiment <- function(FUN, J = 20, tol = 1.0e-06) {
  N_ <- 1000
  x <- runif(N_, 0, 1)
  W <- rep(1, times = N_)
  muBeta0 <- 0
  y <- rbinom(length(x),1,pnorm(FUN(x) + W))
  SigmaBeta0 <- matrix(1,nrow=1,ncol=1)
  fit <- VB_experiment(x, y, W, muBeta0 = muBeta0, SigmaBeta0 = SigmaBeta0, w0 = 1, r0s = 0.01, s0s = 0.01, r0t = 0.01, s0t = 0.01, J = J, tol = 1.0e-06)
  
  varphi <- sqrt(2) * cos(outer(x, pi * (1:J)))
  unknownFunction <- W %*% fit$muBetaq + varphi %*% fit$muThetaq
  ord <- order(x)
  plot(x[ord], unknownFunction[ord], col = 'purple', xlab = 'index', ylab = 'observed/fitted', main = 'Simulation result', ylim = c(-10, 10), type = 'l')
  curve(FUN, from = 0, to = 1, col = 'red', lty = 2, add = TRUE)
  legend('topright', legend = c('fitted values', 'true function'), col = c('purple', 'red'), lty = c(1, 2), bg = 'gray95')
  list(fit = fit, y = y, x = x, varphi = varphi)
}