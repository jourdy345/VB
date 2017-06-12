path=getwd()
source(paste(path,'/GP_code.R',sep=''))

### Data loading
# y
y_tr=scan(paste(path,'/PendulumData/T_tr.txt',sep=''))
ntr=length(y_tr)
y_tst=scan(paste(path,'/PendulumData/T_tst.txt',sep=''))
ntst=length(y_tst)

# X
d=9
X_tr=matrix(scan(paste(path,'/PendulumData/X_tr.txt',sep='')),nr=ntr,nc=d,byrow=T)
X_tst=matrix(scan(paste(path,'/PendulumData/X_tst.txt',sep='')),nr=ntr,nc=d,byrow=T)

# GP_variationial
m=10
S <- rmvnorm(m,mean=rep(0,d),sigma=diag(d))
T=Tfunc(X=X_tr,S=S)$T
fit=VARC(y=y_tr,X=X_tr,T=T,As=25,Ag=25,mul0=rep(0,d),sigl0=10*diag(d),fac=1.5)
names(fit)

# par(mfrow=c(1,2))
plot(fit$lbrecord[,2])
mulq  <- fit$mulq
siglq <- fit$siglq
muaq  <- fit$muaq
EqZ   <- MZ(n,m,T,mulq,siglq)
y_hat <- 