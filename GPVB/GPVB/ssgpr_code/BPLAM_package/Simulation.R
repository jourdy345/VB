library(MASS)

data.gene<-function(n,p)       
  {
    cov<-matrix((0.5)^abs(outer(seq(1,p,1),seq(1,p,1),FUN="-")),ncol=p);
    X<-mvrnorm(n,rep(0,p),cov);X<-pnorm(X);
    fx1<-matrix(sin(2*pi*X[,1])/(2-sin(2*pi*X[,1])),ncol=1); 
    fx2<-matrix(5*X[,2]*(1-X[,2]),ncol=1);
    fx3<-matrix(2*X[,3],ncol=1);
    fx4<-matrix(X[,4],ncol=1); 
    fx5<-matrix(-X[,5],ncol=1);
    epsilon<-matrix(rnorm(n,0,0.5)*(0.5+X[,2]),ncol=1);
    Y<-fx1+fx2+fx3+fx4+fx5+epsilon;
    return(list(X=X,Y=Y));
 }

#---------configeration---------#

   l<-0.2;k <- length(seq(0, 1, l)) - 2
   q <- 3;K <- q + k - 1
   n=200;p=10
   iter=80;sample_size=1000

#---------sample bases-----------#

  set.seed(2014)
  normal_base <- matrix(rnorm (iter*(p+1)*sample_size*2,0,1),nrow=sample_size*2)
  mvnorm_base <- matrix(rnorm (iter*p*sample_size*2*K,0,1),nrow=sample_size*2)
  bernouli_base <- matrix(runif(2*iter*p*sample_size*2),nrow=sample_size*2) 
  gamma_base <- matrix(rgamma (iter*p*sample_size*2,1,1)-1,nrow=sample_size*2)

  subset <- sample(1:(2*sample_size),size=sample_size)
  normb <- normal_base[subset,]
  bernb <- bernouli_base[subset,]
  betab <- mvnorm_base[subset,]
  gammab <- gamma_base[subset,]

#--------------main--------------#

  #source('VB_MC.R')
  # source('VB_Laplace.R')

  data<-data.gene(n,p)
  X<-data$X;Y<-data$Y

  start <- date()
  result<-BPLAM_VB(X,Y,l)         
  end <- date()
  start;end
  #save(result,file='sample_run.result')

#----------variable selection------

   A<-result$p_gamma.a*(1-result$p_gamma.b)
   B<-result$p_gamma.b
   C<-(1-result$p_gamma.a)*(1-result$p_gamma.b)

   A<-t(A);B<-t(B);C<-t(C)
   M<-rbind(B,A,C)

   N1<-sum(apply(M,2,order)[3,]==1)
   N2<-(apply(M,2,order)[3,1]==1)+(apply(M,2,order)[3,2]==1)
   N3<-sum(apply(M,2,order)[3,]==2)
   N4<-(apply(M,2,order)[3,3]==2)+(apply(M,2,order)[3,4]==2)+(apply(M,2,order)[3,5]==2)
   N5<-N1+N3
   N6<-(apply(M,2,order)[3,1]==1||apply(M,2,order)[3,1]==2)+(apply(M,2,order)[3,2]==1||apply(M,2,order)[3,2]==2)+(apply(M,2,order)[3,3]==1||apply(M,2,order)[3,3]==2)+(apply(M,2,order)[3,4]==1||apply(M,2,order)[3,4]==2)+(apply(M,2,order)[3,5]==1||apply(M,2,order)[3,5]==2)
   VS<-c(N1,N2,N3,N4,N5,N6);VS   #should be 2 2 3 3 5 5

