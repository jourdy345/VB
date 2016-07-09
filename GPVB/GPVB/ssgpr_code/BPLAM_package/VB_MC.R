
tpower<-function(x,t)
  {
     q<-3
     a<--(1-t)^(q+1)*6*(q+2*t)/(q+1)/(q+2)
     b<-(1-t)^(q+1)*2*(q+3*t-1)/(q+1)/(q+2)
     (x-t)^q*(x>t)+a*x+b
  }
tpspline<-function(X,l)
  {
     q<-3;
     n<-length(X)
     part1<-matrix(outer(c(X),seq(2,(q-1),1),FUN="^"),ncol=q-2)-X+1/6
     part2<-matrix(outer(c(X),seq(0,(1-l),l),FUN="tpower"),nrow=n);
     B<-matrix(cbind(part1,part2),nrow=n);
     return(B=B);
  }

BPLAM_VB<-function (X,Y,l) 
{
    s <- 1                         #iteration counter
    n <- dim(X)[1];p <- dim(X)[2]
    k <- length(seq(0, 1, l)) - 2
    q <- 3; K<- q + k - 1
    a1 <- 0.5;a2 <- 0.5
    A_delta0.sq <- 1
    B_delta0.sq <- 1

    A_sigma.sq <-rep(1,p);
    B_sigma.sq <-rep(1,p);
    A_tau.sq <-rep(1,p);
    B_tau.sq <-rep(1,p);
    p_gamma.a <-runif(p);
    p_gamma.b <-runif(p);
    mu_u <-runif(1); 
    xi_u <-runif(1);
    mu_a <-runif(p);  
    xi_a <-rep(1,p);
    mu_b <-matrix(runif(K*p),nrow=K,ncol=p);
    xi_b <-matrix(0,nrow=K,ncol=p*K)
    xi_b[matrix(c(rep(1:K,p),(1:(K*p))),ncol=2)]<-1
    B0=X-t(matrix(outer(matrix(apply(X,2,mean),nrow=1),rep(1,n),FUN="*"),ncol=n))
    vec <- c(mu_u,xi_u,mu_a,xi_a,p_gamma.a,mu_b,xi_b,p_gamma.b,A_sigma.sq,
           B_sigma.sq,A_tau.sq,B_tau.sq,A_delta0.sq,B_delta0.sq)
    error=100;indicator=0         #terminate the algorithm if indicator becomes nonzero

    omega <- matrix(NA, ncol = q + k - 1, nrow = k + q - 1)
    for (i in 1:(q - 2)) {
        for (j in 1:(q - 2)) {
            f <- function(x) {
                x^(i - 1) * x^(j - 1)
            }
            omega[i, j] <- i * j * (i + 1) * (j + 1)/(i + j - 
                1) * integrate(f, 0, 1)$value
        }
    }
    for (i in (q - 1):(q + k - 1)) {
        for (j in 1:(q - 2)) {
            f <- function(x) {
                (x - l * (i - q + 1))^(q - 2) * (x^(j - 1))
            }
            omega[i, j] <- (j + 1) * j * q * (q - 1) * integrate(f, 
                l * (i - q + 1), 1)$value
            omega[j, i] <- omega[i, j]
        }
    }
    for (i in (q - 1):(q + k - 1)) {
        for (j in (q - 1):(q + k - 1)) {
            f <- function(x) {
                (x - l * (i - q + 1))^(q - 2) * (x - l * (j - 
                  q + 1))^(q - 2)
            }
            omega[i, j] <- q^2 * (q - 1)^2 * integrate(f, max(l * 
                (i - q + 1), l * (j - q + 1)), 1)$value
        }
    }
    omega<-omega/100

#--------------initialization---------------#

    gamma.a <- outer(rep(1,sample_size),p_gamma.a, FUN='*')
    gamma.a [gamma.a - bernb[,1:p] > 0] <- 1
    gamma.a [gamma.a - bernb[,1:p] < 0] <- 0

    gamma.b <- outer(rep(1,sample_size),p_gamma.b, FUN='*')
    gamma.b [gamma.b - bernb[,(1+iter*p):(p+iter*p)] > 0] <- 1
    gamma.b [gamma.b - bernb[,(1+iter*p):(p+iter*p)] < 0] <- 0

    beta <- betab[,((s-1)*p*K+1):((s-1)*p*K+K)]%*%chol(xi_b[,(1:K)]) +
            outer(rep(1,sample_size),mu_b[,1],FUN='*')
    for (i in 2:p)
     {
       beta <- cbind(beta,betab[,((s-1)*p*K+(i-1)*K+1):((s-1)*p*K+i*K)] %*%
               chol(xi_b[,(K*(i-1)+1):(K*i)]) + outer(rep(1,sample_size),mu_b[,i],FUN='*'))
     }

#------------------y_stars-------------------#

    p1 <- B0 * t(matrix(outer(p_gamma.a*mu_a, rep(1, n), FUN = "*"),ncol = n))
    Ystar_part1 <- matrix(apply(p1, 1, sum), ncol = 1)

    Ystar_part2 <- matrix(rep(0, n), ncol = 1)
    for (i in 1:p) {
        Bj.Betaj <- tpspline(X[, i], l) %*% (mu_b[, i]*p_gamma.b[i])
        Ystar_part2 <- Ystar_part2 + Bj.Betaj
      }

    Ystar <- Y - mu_u - Ystar_part1 - Ystar_part2

    Ytemp_part2 <- matrix(rep(0, n*sample_size), ncol = sample_size)
     for (i in 1:p)
      {
        Bj.Betaj <- tpspline(X[, i], l) %*% t(beta[,(K*(i-1)+1):(K*i)] * 
                     outer(gamma.b[,i],rep(1,K),FUN='*'))
        Ytemp_part2 <- Ytemp_part2 + Bj.Betaj
      }

#------------------- iteration-------------------------#

   while (indicator==0) {

#--------------inverse.gamma's------------------#
 
   delta0.sq <- gammab[,((s-1)*p+sample(1:p,1))]*sqrt(A_delta0.sq/B_delta0.sq^2)+A_delta0.sq/B_delta0.sq

   sigma.sq <- gammab[,((s-1)*p+1):(s*p)]* t(outer(sqrt(A_sigma.sq/B_sigma.sq^2),rep(1,
               sample_size),FUN='*')) + t(outer(A_sigma.sq/B_sigma.sq,rep(1,sample_size),FUN='*'))

   tau.sq <- gammab[,((s-1)*p+1):(s*p)] * t(outer(sqrt(A_tau.sq/B_tau.sq^2),rep(1,sample_size
             ),FUN='*')) + t(outer(A_tau.sq/B_tau.sq,rep(1,sample_size),FUN='*'))

#---------------------mu----------------------#

    Yt <- Ystar + mu_u
    mu_u <- mean(Yt)
    xi_u <- c(B_delta0.sq/(n*A_delta0.sq))
    Ystar <- Yt -mu_u
    mu <- normb[,(iter*p+s)]*sqrt(xi_u)+mu_u

#------------------alpha----------------------#
      
    for (i in 1:p) 
     {
       Yt <- Ystar + p1[, i]
       xi_a[i] <-1/(t(B0[,i])%*%B0[,i]*c(A_delta0.sq/B_delta0.sq) +
                    c(A_sigma.sq[i]/B_sigma.sq[i]))
       mu_a[i] <- xi_a[i]*c(A_delta0.sq/B_delta0.sq)*(B0[,i]%*%Yt)

   #---update the Ystar---#

       Ystar <- Yt - B0[,i]*p_gamma.a[i]*mu_a[i]
     }
    Ystar_part1 <- Y - Ystar - Ystar_part2 - mu_u

#-----------------gamma.a---------------------#

   #---------Ytemp construction--------#

    alpha <- normb[,((s-1)*p+1):(s*p)] * outer (rep(1,sample_size),sqrt(xi_a), 
             FUN='*') + outer (rep(1,sample_size),mu_a, FUN='*')

    Ytemp_part1 <- B0 %*% t(gamma.a*alpha)      #n*sample_size

    Ytemp <- outer(c(Y),rep(1,sample_size),FUN='*') - t(outer(mu,rep(1,n),FUN='*')) -
             Ytemp_part1 - Ytemp_part2

   #------------------------------------#

    asd <- matrix(rep(1, p * p), ncol = p) - diag(c(rep(1, p)))
 
    for (i in 1:p)
     {
       Ytemp_i <- Ytemp + matrix(B0[,i],ncol=1) %*% matrix(gamma.a[,i]*alpha[,i],nrow=1)  

       h1.p1 <- log(delta0.sq * t(B0[,i])%*%B0[,i]+sigma.sq[,i]) 
       h1.p2 <- delta0.sq^2 / (delta0.sq*t(B0[,i])%*%B0[,i] + sigma.sq[,i])   
       h1.p3 <-(B0[,i] %*% Ytemp_i)^2

       h1 <- mean(h1.p1 - h1.p2 * h1.p3)

       h0 <- (p - gamma.a %*% asd[,i]) / (1 + gamma.a %*% asd[,i])
   
       h_gamma.a <-  exp (h1 /2) * exp(mean(log(h0))) *(B_sigma.sq[i] /
                     exp ( digamma ( A_sigma.sq[i] ))) ^ (1/2)
       p_gamma.a[i] <- 1 / (h_gamma.a + 1)

       gamma.a [rep(p_gamma.a[i],sample_size)-bernb[,(s-1)*p+i]>0,i]<-1
       gamma.a [rep(p_gamma.a[i],sample_size)-bernb[,(s-1)*p+i]<0,i]<-0
     
       Ytemp<-Ytemp_i -  matrix(B0[,i],ncol=1) %*% matrix(gamma.a[,i]*alpha[,i],nrow=1)
     }
   #---update the Ystar---#

    Ystar <- Ystar + Ystar_part1
    p1 <- B0 * t(matrix(outer(p_gamma.a*mu_a, rep(1, n), FUN = "*"),ncol = n))
    Ystar_part1 <- matrix(apply(p1, 1, sum), ncol = 1)
    Ystar <- Ystar - Ystar_part1

#------------------beta-----------------#

    for (i in 1:p)
     {
       Bj<-tpspline(X[, i], l)
       Yt <- Ystar +  Bj %*% (mu_b[, i] * p_gamma.b[i])
       xi_b[,(K*(i-1)+1):(K*i)] <- solve(c(A_delta0.sq/B_delta0.sq) * (t(Bj) %*% Bj) + 
                                    A_tau.sq[i]/B_tau.sq[i] * omega)
       mu_b[,i] <- c(A_delta0.sq/B_delta0.sq)*(xi_b[,(K*(i-1)+1):(K*i)] %*% t(Bj) %*%Yt)     

  #---update the Ystar---#

       Ystar <- Yt - Bj %*% (mu_b[, i] * p_gamma.b[i])   
     }

    Ystar_part2 <- Y - Ystar - Ystar_part1 - mu_u

#-------------gamma.b-------------------#

   #-------sample gamma.a---------#

   gamma.a <- outer(rep(1,sample_size),p_gamma.a, FUN='*')
   gamma.a [gamma.a - bernb[,((s-1)*p+1):(s*p)] > 0] <- 1
   gamma.a [gamma.a - bernb[,((s-1)*p+1):(s*p)] < 0] <- 0

   #-------Ytemp construction------#

   Ytemp_part1 <- B0 %*% t(gamma.a*alpha)      #n*sample_size
     
   beta <- betab[,((s-1)*p*K+1):((s-1)*p*K+K)]%*%chol(xi_b[,(1:K)]) +
            outer(rep(1,sample_size),mu_b[,1],FUN='*')
    for (i in 2:p)
     {
       beta <- cbind(beta,betab[,((s-1)*p*K+(i-1)*K+1):((s-1)*p*K+i*K)] %*%
               chol(xi_b[,(K*(i-1)+1):(K*i)]) + outer(rep(1,sample_size),mu_b[,i],FUN='*'))
     }
    
   Ytemp_part2 <- matrix(rep(0, n*sample_size), ncol = sample_size)
     for (i in 1:p)
      {
       Bj.Betaj <- tpspline(X[, i], l) %*% t(beta[,(K*(i-1)+1):(K*i)] * 
                   outer(gamma.b[,i],rep(1,K),FUN='*'))
       Ytemp_part2 <- Ytemp_part2 + Bj.Betaj
      }
   Ytemp <- outer(c(Y),rep(1,sample_size),FUN='*') - t(outer(mu,rep(1,n),FUN='*')) -
             Ytemp_part1 - Ytemp_part2

  #-------------------------------#

   asd <- matrix(rep(1, p * p), ncol = p) - diag(c(rep(1, p)))

   for (i in 1:p) 
    {
      Bj<-tpspline(X[, i], l)
      Ytemp_i<-Ytemp +  Bj%*%t(beta[,(K*(i-1)+1):(K*i)]*outer(gamma.b[,i],rep(1,K),FUN='*'))

      h1.p1 <- rep(0,sample_size) ; h1.p2 <- rep(0,sample_size)
        for (t in 1:sample_size)
         {
           h1.p1[t] <- determinant(delta0.sq[t]*t(Bj)%*%Bj + tau.sq[t,i] * 
                    omega, logarithm = TRUE)$modulus
           h1.p2[t] <- delta0.sq[t] * Ytemp_i[,t] %*% Bj %*% solve(delta0.sq[t]*t(Bj) %*% 
                    Bj + tau.sq[t,i] * omega) %*% t(Bj) %*% (delta0.sq[t]*Ytemp_i[,t])
         } 
      h1 <- mean(h1.p1-h1.p2)
      h0 <- (p - gamma.b %*% asd[,i]) / (1 + gamma.b %*% asd[,i])    
      h_gamma.b <- (B_tau.sq[i]/(K*determinant(omega, logarithm = F)$modulus *
                      exp(digamma(A_tau.sq[i]))))^(1/2) * exp(h1 /2)*exp(mean(log(h0)))
      p_gamma.b[i] <- 1 / (h_gamma.b + 1)
      gamma.b [rep(p_gamma.b[i],sample_size)-bernb[,(s-1)*p+i+iter*p]>0,i]<-1
      gamma.b [rep(p_gamma.b[i],sample_size)-bernb[,(s-1)*p+i+iter*p]<0,i]<-0
      Ytemp<-Ytemp_i -  Bj%*%t(beta[,(K*(i-1)+1):(K*i)]*outer(gamma.b[,i],rep(1,K),FUN='*'))
    }     

  #---update the Ystar---#

    Ystar <- Ystar + Ystar_part2
    Ystar_part2 <- matrix(rep(0, n), ncol = 1)
    for (i in 1:p) 
     {
       Bj.Betaj <- tpspline(X[, i], l) %*% (mu_b[, i]*p_gamma.b[i])
       Ystar_part2 <- Ystar_part2 + Bj.Betaj
      }
    Ystar <- Ystar - Ystar_part2

#---------------sigma.sq--------------------#

    A_sigma.sq <- a1 + p_gamma.a / 2
    B_sigma.sq <- a2 + p_gamma.a * (xi_a + mu_a^2) / 2

#-----------------tau.sq--------------------#

    A_tau.sq <- a1 + p_gamma.b * K / 2
      for (i in 1:p)
       {
         B_tau.sq[i] <- a2 + p_gamma.b[i] * sum(diag(t(omega)%*%(mu_b[,i]%*%t(mu_b[,i]) +
                        xi_b[,(K*(i-1)+1):(K*i)]))) / 2
       }

#----------------delta0.sq-------------------#

    A_delta0.sq <- a1 + n / 2

  #----------sample gamma.b------------#

    gamma.b <- outer(rep(1,sample_size),p_gamma.b, FUN='*')
    gamma.b [gamma.b - bernb[,((s-1)*p+1+iter*p):(s*p+iter*p)] > 0] <- 1
    gamma.b [gamma.b - bernb[,((s-1)*p+1+iter*p):(s*p+iter*p)] < 0] <- 0

    Ytemp_part2 <- matrix(rep(0, n*sample_size), ncol = sample_size)
      for (i in 1:p)
        {
         Bj.Betaj <- tpspline(X[, i], l) %*% t(beta[,(K*i-5):(K*i)] * 
                     outer(gamma.b[,i],rep(1,K),FUN='*'))
         Ytemp_part2 <- Ytemp_part2 + Bj.Betaj
        }
    Ytemp <- outer(c(Y),rep(1,sample_size),FUN='*') - t(outer(mu,rep(1,n),FUN='*')) -
             Ytemp_part1 - Ytemp_part2

    B_delta0.sq <- a2 + mean(diag(t(Ytemp)%*%Ytemp)) / 2
 
    vec_temp <- c(mu_u,xi_u,mu_a,xi_a,p_gamma.a,mu_b,xi_b,p_gamma.b,A_sigma.sq,
                B_sigma.sq,A_tau.sq,B_tau.sq,A_delta0.sq,B_delta0.sq)
    error <- sqrt(mean((vec_temp/vec-1)^2))
    vec <- vec_temp
    cat(s, "\n")
    s=s+1;indicator=(error<10^(-3))+(s>80)
  }
return(list(mu_u=mu_u,xi_u=xi_u,mu_a=mu_a,xi_a=xi_a,mu_b=mu_b,xi_b=xi_b,p_gamma.a=p_gamma.a,
p_gamma.b=p_gamma.b,A_sigma.sq=A_sigma.sq,B_sigma.sq=B_sigma.sq,A_tau.sq=A_tau.sq,
B_tau.sq=B_tau.sq,A_delta.sq=A_delta0.sq,B_delta.sq=B_delta0.sq))
}

