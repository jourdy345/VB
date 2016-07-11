require('gpr') # for 'minimize'
X_tr <- read.table(file = '~/Desktop/Github/VB/GPVB/GPVB/PendulumData/X_tr.txt')
X_tst <- read.table(file = '~/Desktop/Github/VB/GPVB/GPVB/PendulumData/X_tst.txt')
T_tr <- read.table(file = '~/Desktop/Github/VB/GPVB/GPVB/PendulumData/T_tr.txt')
T_tst <- read.table(file = '~/Desktop/Github/VB/GPVB/GPVB/PendulumData/T_tst.txt')
loghyp <- rep(1, 11)
X_tr <- unname(as.matrix(X_tr))
X_tst <- unname(as.matrix(X_tst))
T_tr <- unname(as.matrix(T_tr))
T_tst <- unname(as.matrix(T_tst))


source('~/Desktop/Github/VB/cosine/demo/vbgpspectral.R')

minim <- function (X, f, .length, x, y)
{
    toprint = FALSE
    INT = 0.1
    EXT = 3
    MAX = 20
    RATIO = 10
    SIG = 0.1
    RHO = SIG/2
    if (is.array(.length)) {
        if (max(dim(.length)) == 2) {
            red = .length[2]
            .length = .length[1]
        }
    }
    else {
        red = 1
    }
    if (.length > 0) {
        S = "Linesearch"
    }
    else {
        S = "Function evaluation"
    }
    i.m = 0
    ls_failed = 0
    f_out = eval(call(f, X, x, y))
    f0 = c(f_out[1][[1]])
    df0 = c(f_out[2][[1]])
    fX = f0
    i.m = i.m + (.length < 0)
    s = -df0
    s = round(s * 100000)/100000
    d0 = -c(crossprod(s))
    x3 = red/(1 - d0)
    mainloop = TRUE
    while (i.m < abs(.length) && mainloop) {
        i.m = i.m + (.length > 0)
        X0 = X
        F0 = f0
        dF0 = df0
        if (.length > 0) {
            M = MAX
        }
        else {
            M = min(MAX, -.length - i.m)
        }
        whilerun = TRUE
        while (whilerun == TRUE) {
            x2 = 0
            f2 = f0
            d2 = d0
            f3 = f0
            df3 = df0
            success = FALSE
            while (success == FALSE && M > 0) {
                M = M - 1
                i.m = i.m + (.length < 0)
                options(show.error.messages = FALSE)
                f_out2 = eval(call(f, c(X + x3[1] * s), x, y))
                f3 = c(f_out2[1][[1]])
                df3 = c(f_out2[2][[1]])
                f3 = round(f3 * 100000)/100000
                df3 = round(df3 * 100000)/100000
                if (is.na(f3) || is.infinite(f3) || is.nan(f3) ||
                  any(is.nan(df3) || is.na(df3) || is.infinite(df3))) {
                  cat(" ")
                  x3 = (x2 + x3)/2
                }
                else {
                  success = TRUE
                }
                options(show.error.messages = TRUE)
            }
            if (f3 < F0) {
                X0 = X + x3[1] * s
                F0 = f3
                dF0 = df3
            }
            d3 = c(crossprod(df3, s))
            if (d3 > SIG * d0 || f3 > f0 + x3 * RHO * d0 || M ==
                0) {
                whilerun = FALSE
                break
            }
            x1 = x2
            f1 = f2
            d1 = d2
            x2 = x3
            f2 = f3
            d2 = d3
            A = 6 * (f1 - f2) + 3 * (d2 + d1) * (x2 - x1)
            B = 3 * (f2 - f1) - (2 * d1 + d2) * (x2 - x1)
            x3 = x1 - d1 * (x2 - x1)^2/(B + sqrt(abs(B * B -
                A * d1 * (x2 - x1))))
            if ((B * B - A * d1 * (x2 - x1) < 0)[1] || is.nan(x3) ||
                is.infinite(x3) || x3 < 0) {
                x3 = x2 * EXT
            }
            else if (x3 > x2 * EXT) {
                x3 = x2 * EXT
            }
            else if (x3 < x2 + INT * (x2 - x1)) {
                x3 = x2 + INT * (x2 - x1)
            }
            x3 = round(x3 * 100000)/100000
        }
        while ((abs(d3) > -SIG * d0 || f3 > f0 + x3 * RHO * d0) &&
            M > 0) {
            if (d3 > 0 || f3 > f0 + x3 * RHO * d0) {
                x4 = x3
                f4 = f3
                d4 = d3
            }
            else {
                x2 = x3
                f2 = f3
                d2 = d3
            }
            if (f4 > f0) {
                x3 = x2 - (0.5 * d2 * (x4 - x2)^2)/(f4 - f2 -
                  d2 * (x4 - x2))
            }
            else {
                A = 6 * (f2 - f4)/(x4 - x2) + 3 * (d4 + d2)
                B = 3 * (f4 - f2) - (2 * d2 + d4) * (x4 - x2)
                x3 = x2 + (sqrt(B * B - A * d2 * (x4 - x2)^2) -
                  B)/A
            }
            if (is.nan(x3) || is.infinite(x3)) {
                x3 = (x2 + x4)/2
            }
            x3 = max(min(x3, x4 - INT * (x4 - x2)), x2 + INT *
                (x4 - x2))
            f_out3 = eval(call(f, c(X + x3 * s), x, y))
            f3 = c(f_out3[1][[1]])
            df3 = c(f_out3[2][[1]])
            if (f3 < F0) {
                x3 = x3[[1]]
                X0 = X + x3 * s
                F0 = f3
                dF0 = df3
            }
            M = M - 1
            i.m = i.m + (.length < 0)
            d3 = c(crossprod(df3, s))
        }
        if (abs(d3) < -SIG * d0 && f3 < f0 + x3 * RHO * d0) {
            x3 = x3[[1]]
            X = X + x3 * s
            f0 = f3
            fX = c(fX, f0)
            cat(S, i.m, "; Value ", f0, '\n')
            s = (((c(crossprod(df3)) - c(crossprod(df0, df3)))[1])/((c(crossprod(df0)))[[1]]) * s) - df3
            df0 = df3
            d3 = d0
            d0 = c(crossprod(df0, s))
            if (d0 > 0) {
                s = -df0
                d0 = -c(crossprod(s))
            }
            x3 = x3 * min(RATIO, d3/(d0 - (2^(-1022))))
            ls_failed = 0
        }
        else {
            X = X0
            f0 = F0
            df0 = dF0
            if (ls_failed || i.m > abs(.length)) {
                mainloop = 0
                break
            }
            s = -df0
            d0 = -c(crossprod(s))
            x3 = 1/(1 - d0)
            ls_failed = 1
        }
    }
    return(list(X=X, fX=fX, i.m=i.m))
}


ssgpr <- function(optimizeparams, x_tr, y_tr, x_tst = NULL) {
  n <- dim(x_tr)[1]
  D <- dim(x_tr)[2]
  m <- (length(optimizeparams) - D - 2) / D
  ell <- exp(optimizeparams[1:D])
  sf2 <- exp(2 * optimizeparams[D+1])
  sn2 <- exp(2 * optimizeparams[D+2])
  w <- matrix(c(optimizeparams[(D+3):length(optimizeparams)]), nrow = m, ncol = D)
  # cat('dim of w: ', dim(w), '\n')
  # cat('diag(1/ell): ', diag(ell), '\n')
  # cat('ell: ', ell, '\n')
  if (length(ell) == 1) {
    w <- w / ell # dividing each row of w by ell element-wise... needs to be revised if dimension of w changes
  } else {
    w <- w %*% diag(1 / ell) # divide each row of w by ell element-wise
  }

  phi <- tcrossprod(x_tr, w)
  phi <- cbind(cos(phi), sin(phi))
  
  R <- chol((sf2/m) * crossprod(phi) + sn2 * diag(2 * m))
  PhiRi <- phi %*% solve(R)
  RtiPhit <- t(PhiRi)
  Rtiphity <- RtiPhit %*% y_tr
  
  if (is.null(x_tst)) {
    out1 <- 0.5 / sn2 * (sum(y_tr^2) - sf2 / m * sum(Rtiphity^2)) + sum(log(diag(R))) + (n / 2 - m) * log(sn2) + n / 2 + log(2 * pi)
    
    out2 <- rep(0, D+2+D*m)
    
    A <- cbind(y_tr / sn2 - PhiRi %*% ((sf2 / sn2 / m) * Rtiphity), sqrt(sf2 / sn2 / m) * PhiRi)
    diagfact <- -1 / sn2 + rowSums(A^2)
    Aphi <- crossprod(A, phi)
    B <- ((A %*% Aphi[,1:m]) * phi[, (m+1):ncol(phi)]) - ((A %*% Aphi[, (m+1):ncol(Aphi)]) * phi[,1:m])
    
    for (d in 1:D) {
      out2[d] <- -0.5 * 2 * sf2 / m * (crossprod(x_tr[,d], B) %*% w[,d])
    }
    out2[D+1] <- 0.5 * 2 * (sf2 / m) * (n * m / sn2 - sum(Aphi^2))
    out2[D+2] <- -0.5 * sum(diagfact) * 2 * sn2
    
    for (d in 1:D) {
      out2[(D+2+(d-1)*m+1):(D+2+d*m)] <- 0.5 * 2 * sf2 / m * (c(crossprod(x_tr[,d], B)) / ell[d])
    }
  } else {
    ns <- dim(x_tst)[1]
    out1 <- out2 <- rep(0, ns)
    alfa <- sf2 / m * (solve(R, Rtiphity))
    
    chunksize <- 5000
    allxstar <- x_tst
    
    for (beg_chunk in seq(from = 1, to = ns, by = chunksize)) {
      end_chunk <- min(beg_chunk + chunksize - 1, ns)
      x_tst <- allxstar[beg_chunk:end_chunk,]
      
      phistar <- tcrossprod(x_tst, w)
      phistar <- cbind(cos(phistar), sin(phistar))
      out1[beg_chunk:end_chunk] <- phistar %*% alfa
      
      out2[beg_chunk:end_chunk] <- sn2 * (1 + sf2 / m * rowSums((phistar %*% solve(R))^2))
    }
    
  }
  list(out1 = out1, out2 = out2)
}


ssgpr_ui <- function(x_tr, y_tr, x_tst, y_tst, m, iteropt = NULL, loghyper = NULL) {
  meanp <- mean(apply(y_tr, 2, mean))
  y_tr <- scale(y_tr, scale = FALSE)
  n <- dim(x_tr)[1]
  D <- dim(x_tr)[2]
  if ((!is.null(loghyper)) & (length(loghyper) == (D+2))) {
    lengthscales <- loghyper[1:D]
    covpower <- loghyper[D+1]
    noisepower <- loghyper[D+2]
    nlml <- Inf
    optimizeparams <- c(lengthscales, covpower, noisepower)
    for (k in 1:100) {
      otherparams <- rnorm(m * D)
      nlmlc <- ssgpr(c(optimizeparams, otherparams), x_tr, y_tr)$out1
      if (nlmlc < nlml) {
        w_save <- otherparams
        nlml <- nlmlc
      }
    }
    otherparams <- w_save
  } else if (is.null(loghyper)) {
    lengthscales <- log((apply(x_tr, 2, max) - apply(x_tr, 2, min)) / 2)
    lengthscales[lengthscales < -1e2] <- -1e2
    covpower <- 0.5 * log(apply(y_tr, 2, function(x) mean((x-mean(x))^2)))
    noisepower <- 0.5 * log(apply(y_tr, 2, function(x) mean((x-mean(x))^2)) / 4)
    nlml <- Inf
    optimizeparams <- c(lengthscales, covpower, noisepower)
    for (k in 1:100) {
      otherparams <- rnorm(m * D)
      nlmlc <- ssgpr(c(optimizeparams, otherparams), x_tr, y_tr)$out1
      if (nlmlc < nlml) {
        w_save <- otherparams
        nlml <- nlmlc
      }
    }
    otherparams <- w_save
  } else if ((!is.null(loghyper)) & (length(loghyper) != D+2+D*m)) {
    stop('Incorrect number of hyperparameters.')
  } else if ((!is.null(loghyper)) & (length(loghyper) == D+2+D*m)) {
    optimizeparams <- loghyper[1:(D+2)]
    otherparams <- loghyper[(D+3):(D+2+D*m)]
  }
  if (is.null(iteropt)) iteropt <- -1000
  
  optimizeparams <- c(optimizeparams, otherparams)
  temp <- minim(optimizeparams, 'ssgpr', iteropt, x_tr, y_tr)

  optimizeparams <- temp$X
  convergence <- temp$fX
  loghyper <- optimizeparams

  res <- ssgpr(optimizeparams, x_tr, y_tr, x_tst)

  mu <- res$out1
  S2 <- res$out2
  mu <- mu + meanp  
  NMSE <- mean((mu-y_tst)^2) / mean((meanp - y_tst)^2)
  
  NMLP <- -0.5 * mean((-(mu - y_tst)^2) / S2 - log(2 * pi) - log(S2))
  list(NMSE = NMSE, mu = mu, S2 = S2, NMLP = NMLP, loghyper = loghyper, convergence = convergence)
}

# # fit <- ssgpr_ui(X_tr, T_tr, X_tst, T_tst, 100, -1000, loghyp)
# fit <- ssgpr_ui(X_tr, T_tr, X_tr, T_tr, 100, -1000, loghyp)

# library(glmnet)
# cvfit <- glmnet::cv.glmnet(X_tr, T_tr)
# coef(cvfit, s = 'lambda.1se')

# # Set up prior parameters
# J <- 20
# rsig.0<-0.01
# ssig.0<-0.01
# rtau.0<-0.01
# stau.0<-0.01
# w0<-1
# #mubeta.0<-c(0,0)
# #sigbeta.0<-diag(2)
# mubeta.0<-0
# sigbeta.0<-matrix(1,nrow=1,ncol=1)
# prior.parms<-list(rsig.0=rsig.0,ssig.0=ssig.0,rtau.0=rtau.0,stau.0=stau.0,w0=w0,mubeta.0=mubeta.0,sigbeta.0=sigbeta.0)
# fit2 <- vbgpspectral(T_tr, X_tr[,7], rep(1, length(T_tr)), 20, 1.0e-05, prior.parms = prior.parms, mupsi.q.start = 1)
# x <- X_tr[,7]
# vphi<-sqrt(2)*cos(outer(x,pi*(1:J)))
# fitted2<-(vphi[,1:length(fit2$mutheta.q)]%*%fit2$mutheta.q)
# # fitted2<-fitted2-mean(fitted2)
# # fitted2 <- vphi %*% fit2$mutheta.q + fit2$mubeta.q
# o <- order(X_tr[,7])
# plot(x[o], T_tr[o], type = 'l', lwd = 2, lty = 6, col = 'darkgreen')
# lines(x[o], fitted2[o], lwd = 2, lty = 3, col = 'red')
# lines(x[o], fit$mu[o], lwd = 2, lty = 1)



# set.seed(1)
# n<-500
# J <- 20
# x<-0.1+0.8*runif(n)
# vphi2 <- sqrt(2)*cos(outer(x,pi*(1:J)))
# loghyper = rep(1,3)
# #Z<-cbind(rep(1,times=n),runif(n))
# #y<-sin(x*pi)+Z%*%c(1,1)+0.1*rnorm(n)
# Z<-rep(1,times=n)
# yStar <- sin(x*pi)+Z
# y<-yStar+0.1*rnorm(n)
# rsig.0<-0.01
# ssig.0<-0.01
# rtau.0<-0.01
# stau.0<-0.01
# w0<-1
# #mubeta.0<-c(0,0)
# #sigbeta.0<-diag(2)
# mubeta.0<-0
# sigbeta.0<-matrix(1,nrow=1,ncol=1)
# prior.parms<-list(rsig.0=rsig.0,ssig.0=ssig.0,rtau.0=rtau.0,stau.0=stau.0,w0=w0,mubeta.0=mubeta.0,sigbeta.0=sigbeta.0)

# fit_SSGP <- c()
# for (i in seq(from = 20, to = 200, by = 10)) {
#   fit_SSGP <- c(fit_SSGP, (ssgpr_ui(as.matrix(x), as.matrix(y), as.matrix(x), as.matrix(yStar), i, -1000, rep(1, 3)))$NMSE)
# }

# fit_BSAR <- c()
# for (j in seq(from = 20, to = 200, by = 10)) {
#   temp <- vbgpspectral(y, x, Z = rep(1, n), T = i, tol = 1.0e-05, prior.parms = prior.parms, mupsi.q.start = 1)
#   fitted_BSAR <- vphi2[,1:length(temp$mutheta.q)]%*%temp$mutheta.q
#   res <- mean((fitted_BSAR - yStar)^2) / mean((yStar - mean(y))^2)
#   fit_BSAR <- c(fit_BSAR, res)
# }

# temp <- vbgpspectral(y, x, Z = rep(1, n), T = 20, tol = 1.0e-05, prior.parms = prior.parms, mupsi.q.start = 1)
# fitted_BSAR <- vphi2[,1:length(temp$mutheta.q)]%*%temp$mutheta.q
# res <- mean((fitted_BSAR - yStar)^2) / mean((yStar - mean(y))^2)

# plot(seq_along(fit_SSGP), fit_SSGP, pch = 1, ylim = range(c(fit_SSGP, fit_BSAR)), col = 'darkgreen', xlab = '# of basis functions', ylab = 'NMSE', main = 'SSGP vs BSAR')
# lines(seq_along(fit_SSGP), fit_SSGP, lty = 2, col = 'darkgreen')
# points(seq_along(fit_BSAR), fit_BSAR, pch = 2, col = 'purple')
# lines(seq_along(fit_BSAR), fit_BSAR, lty = 3, col = 'purple')


# t1 <- Sys.time()
# fitt <- ssgpr_ui(as.matrix(x), as.matrix(y), as.matrix(x), as.matrix(y), 100, -100, rep(1, 3))
# t2 <- Sys.time()
# z1 <- difftime(t2, t1)
# class(z1) <- NA
# z1 <- round(z1[1], digits = 4)
# t3 <- Sys.time()
# fitt2 <- vbgpspectral(y, x, Z = rep(1, n), T = 20, tol = 1.0e-05, prior.parms = prior.parms, mupsi.q.start = 1)
# t4 <- Sys.time()
# z2 <- difftime(t4, t2)
# class(z2) <- NA
# z2 <- round(z2[1], digits = 4)
# vphi2 <- sqrt(2)*cos(outer(x,pi*(1:J)))
# fitted22 <- vphi2[,1:length(fitt2$mutheta.q)]%*%fitt2$mutheta.q
# fitted22 <- fitted22 - mean(fitted22)
# o <- order(x)
# y2 <- y - mean(y)
# fitmu <- fitt$mu - mean(fitt$mu)
# plot(x[o], y2[o], ylim = range(c(fitted22, fitmu, y2)), xlab = '', ylab = '', main = expression(paste('sin(',pi, 'x)+Z+0.1N(0,1)')))
# lines(x[o], fitted22[o], lwd = 2, lty = 3, col = 'red')
# lines(x[o], fitmu[o], lwd = 2, lty = 6, col = 'darkgreen')
# legend('topright', lty = c(NA, 3, 6), pch = c(1, NA, NA),  col = c(1, 'red', 'darkgreen'), legend = c('true', paste('BSAR =',z2,'s'), paste('SSGP=',z1,'s')), bg = 'gray90')



compareSSGPvsBSAR <- function(data = 'pendulum', fit = 'training') {
  ############################################################################################
  ############################################################################################
  ############     data: which data,                                              ############
  ############       c('pendulum', 'elevators', 'kin', 'pol', 'pumadyn', 'simul') ############
  ############                                                                    ############
  ############     fit: return training MSE or test MSE?                          ############
  ############       c('training', 'test')                                        ############
  ############################################################################################
  ############################################################################################
  switch()
}
#########################################################
##########                                ###############
##########  Pendulum data: Training data  ###############
##########                                ############### 
#########################################################
y <- T_tr
x <- X_tr[,7]
x <- pnorm(x)
x_rest <- X_tr[,-7]
J <- 20
rsig.0<-0.01
ssig.0<-0.01
rtau.0<-0.01
stau.0<-0.01
w0<-1
mubeta.0<-rep(0, times = ncol(x_rest))
sigbeta.0 <- diag(1, length(mubeta.0))
prior.parms<-list(rsig.0=rsig.0,ssig.0=ssig.0,rtau.0=rtau.0,stau.0=stau.0,w0=w0,mubeta.0=mubeta.0,sigbeta.0=sigbeta.0)
t1 <- Sys.time()
fitt <- ssgpr_ui(X_tr, as.matrix(y), X_tr, as.matrix(y), 100, -100, rep(1, 11))
t2 <- Sys.time()
z1 <- difftime(t2, t1)
class(z1) <- NA
z1 <- round(z1[1], digits = 4)
t3 <- Sys.time()
fitt2 <- vbgpspectral(y, x, Z = x_rest, T = 20, tol = 1.0e-05, prior.parms = prior.parms, mupsi.q.start = 1)
t4 <- Sys.time()
z2 <- difftime(t4, t3)
class(z2) <- NA
z2 <- round(z2[1], digits = 4)
vphi2 <- sqrt(2) * cos(outer(x, pi * (1:J)))
fitted22 <- vphi2[,1:length(fitt2$mutheta.q)] %*% fitt2$mutheta.q
fitted22 <- fitted22 - mean(fitted22)
o <- order(x)
y2 <- y - mean(y)
fitmu <- fitt$mu - mean(fitt$mu)
plot(x[o], y2[o], ylim = range(c(fitted22, fitmu, y2)), xlab = '', ylab = '', main = 'Pendulum training data')
lines(x[o], fitted22[o], lwd = 2, lty = 3, col = 'red')
lines(x[o], fitmu[o], lwd = 2, lty = 6, col = 'darkgreen')
legend('topright', lty = c(NA, 3, 6), pch = c(1, NA, NA),  col = c(1, 'red', 'darkgreen'), legend = c('true', paste('BSAR =',z2,'s'), paste('SSGP=',z1,'s')), bg = 'gray90')




# LASSOselectedVariable <- glmnet::cv.glmnet(elevator_X_tr, c(elevator_T_tr))
# coefficients <- coef(LASSOselectedVariable, s = "lambda.1se")
