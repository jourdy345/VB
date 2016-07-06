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
    print('min>')
    f_out = eval(call(f, X, x, y))
    f0 = c(f_out[1][[1]])
    df0 = c(f_out[2][[1]])
    fX = f0
    i.m = i.m + (.length < 0)
    s = -df0
    s = round(s * 10000)/10000
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
                print('min>>')
                f_out2 = eval(call(f, c(X + x3[1] * s), x, y))
                f3 = c(f_out2[1][[1]])
                df3 = c(f_out2[2][[1]])
                f3 = round(f3 * 10000)/10000
                df3 = round(df3 * 10000)/10000
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
            x3 = round(x3 * 10000)/10000
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
            print('min>>>')
            f_out3 = eval(call(f, c(X + x3 * s), x, y))
            f3 = f_out3[1][[1]][[1]]
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
            fX = t(cbind(t(fX), f0))
            s = (((c(crossprod(df3)) - c(crossprod(df0, df3)))[1])/((c(crossprod(df0)))[[1]]) * s) - df3
            df0 = df3
            d3 = d0
            d0 = t(df0) %*% s
            if (d0 > 0) {
                s = -df0
                d0 = -t(s) %*% s
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
  w <- w %*% diag(1 / ell) # divide each row of w by ell element-wise
  
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
    ns <- dim(x_tst)[2]
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
  meanp <- apply(y_tr, 2, mean)
  y_tr <- scale(y_tr, scale = FALSE)
  n <- dim(x_tr)[1]
  D <- dim(x_tr)[2]
  
  if ((!is.null(loghyper)) & (length(loghyper) == D+2)) {
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
  print('>')
  if (is.null(iteropt)) iteropt <- -1000
  
  print('>>')
  optimizeparams <- c(optimizeparams, otherparams)
  temp <- minim(optimizeparams, 'ssgpr', iteropt, x_tr, y_tr)
  print('>>>')
  optimizeparams <- temp$X
  convergence <- temp$fX
  loghyper <- optimizeparams
  print('>>>>')
  res <- ssgpr(optimizeparams, x_tr, y_tr, x_tst)
  print('>>>>>')
  mu <- res$out1
  S2 <- res$out2
  
  NMSE <- mean((mu-y_tst)^2) / mean((meanp - y_tst)^2)
  
  NMLP <- -0.5 * mean((-(mu - y_tst)^2) / S2 - log(2 * pi) - log(S2))
  list(NMSE = NMSE, mu = mu, S2 = S2, NMLP = NMLP, loghyper = loghyper, convergence = convergence)
}

fit <- ssgpr_ui(X_tr, T_tr, X_tst, T_tst, 100, -1000, loghyp)


