require(truncnorm)
require(mvtnorm)
require(HyperbolicDist)
require(crch)
require(MCMCpack)

# /*
# ****************************************************************
# * RNDGAMLO.G
# *
# *  f(x) \propto I( x> low) x^{alpha-1}exp(-beta*x)
# *
# *	INPUT:
# *		nr		= number of rows
# *		nc		= number of columns
# *		alpha	= shape parameters of gamma: scalar 
# *		beta		= scale parameter: scalar
# *		low		= lower limit: scalar.
# *
# *	OUTPUT
# *		truncated gamma
# ******************************************************************
# */
'rndgamlo'=function(alpha,beta,low)
{
  a2	= alpha
  glow = pgamma(beta*low,shape=alpha);
  if (glow >= 0.9999){   #@ lower bound is big: use shifted exponential approximation in tail @
    x = low-log(1-runif(1))/beta;
  }else{
    p=runif(1)*(1-glow)+glow;
    x=qgamma(p,shape=a2)/beta;	
  }
  x=(low-x)*(x<low)+x;
  return(x);
}


#Transform lower triangle matrix (including diagonal) to a column vector
'vech'=function(mat){
  nr=dim(mat)[1]
  nc=dim(mat)[2]
  vec=NULL
  for(i in 1:min(nr,nc)){
    veci=mat[i,1:i]
    vec=c(vec,veci)
  }
  vec
}



# @ Computes log IG density @
'LogfIG' = function(x,r0,s0) (r0/2)*log(s0/2)-lgamma(r0/2)-  (r0/2+1)*log(x)-s0/(2*x)

# @ QuadMult:  gets quadratic form x'Q*x where qvech = vech(Q) is given @

'QuadMult' = function(x,qvech,quadfacts)
{
  colSums(quadfacts[,1]*x[quadfacts[,2]]*qvech*x[quadfacts[,3]])
}

######################
# Calcuate Functions #
######################

'GetFreef' = function(theta,phixobs,phixgrid,nbasis, nobs, ngrid)
{	
  fxobs = phixobs%*%theta
  fxgrid = phixgrid%*%theta
  return(list(fxobs=fxobs,fxgrid=fxgrid))
}


'GetUpf' = function(fpm,theta,phixobs,phixgrid,quadfacts,nbasis,nr,nobs,ngrid)
{
  # @ Quadratic Form Multiplication @
  fxobs = QuadMult(theta,phixobs,quadfacts)
  fxgrid = QuadMult(theta,phixgrid,quadfacts)
  
  # @ Make it plus or minus @
  fxobs = fpm*fxobs
  fxgrid = fpm*fxgrid
  
  return(list(fxobs=fxobs, fxgrid=fxgrid))
}


'GetConvexf' = function(fpm,alpha,theta,xobs,xgrid,xmid,phixobs,phixgrid,quadfacts,
						nbasis,nr,nobs,ngrid)
{	
  # @ Quadratic Form Multiplication @
  fxobs = QuadMult(theta,phixobs,quadfacts)
  fxgrid = QuadMult(theta,phixgrid,quadfacts)			
  
  # @ Add linear term @
  fxobs = fpm*fxobs+alpha*(xobs-xmid)
  fxgrid = fpm*fxgrid+alpha*(xgrid-xmid) 
  
  return(list(fxobs=fxobs, fxgrid=fxgrid))
}


'GetConcavef' = function(fpm,alpha,theta,xobs,xgrid,xmid,phixobs,phixgrid,quadfacts,
						nbasis,nr,nobs,ngrid)
{	
  # @ Quadratic Form Multiplication @
  fxobs = QuadMult(theta,phixobs,quadfacts)
  fxgrid = QuadMult(theta,phixgrid,quadfacts)			
  
  # @ make them into increasing, concave functions @
  fxobs=-fxobs
  fxgrid=-fxgrid
  
  # @ Add linear term @
  fxobs = fpm*fxobs+alpha*(xobs-xmid) 
  fxgrid = fpm*fxgrid+alpha*(xgrid-xmid) 
  
  return(list(fxobs=fxobs, fxgrid=fxgrid))
}


'GetSf' = function(fpm,omega,psi,alpha,theta,xobs,xgrid,phixobs,phixgrid,xdelta,xinxgrid,
                 	xidelta,xrange,xmid,quadfacts,intsimpfacts,nbasis,nr,nobs,ngrid)
{
  
  # @ Compute z2 at grid and observations @
  # @ Quadratic Form Multiplication @
  z2xobs = QuadMult(theta,phixobs,quadfacts)
  z2xgrid = QuadMult(theta,phixgrid,quadfacts)
  
  # @ Compute squish function at grid and observations @
  hxobs = SquishDown(xobs,psi,omega,nobs)
  hxgrid = SquishDown(xgrid,psi,omega,ngrid)
  
  # @ Second derivatives at xgrid and xobs @
  f2xobs = z2xobs*hxobs
  f2xgrid = z2xgrid*hxgrid
  
  # @ First derivatives at xgrid and xobs  @
  f1xgrid = InTrapCum(f2xgrid,xdelta,ngrid)
  f1xobs = IntFobs(f2xobs,f2xgrid,f1xgrid,xinxgrid,xidelta,nobs,ngrid)
  
  # @ First derivative needs to be positive. @
  f1min = min(c(0,min(f1xgrid))) 
  
  # @ Compute f at xgrid and xobs @
  fxgrid = InTrapCum(f1xgrid,xdelta,ngrid)
  fxobs = IntFobs(f1xobs,f1xgrid,fxgrid,xinxgrid,xidelta,nobs,ngrid)
  
  # @ Mean center f @
  fint = IntSimpsonfxgrid(fxgrid,xdelta,intsimpfacts,ngrid)
  fxobs = fxobs-fint/xrange	
  fxgrid = fxgrid-fint/xrange
  
  # @ Add linear term @
  fxobs = fpm*(fxobs)+(alpha-f1min)*(xobs-xmid)
  fxgrid = fpm*(fxgrid)+(alpha-f1min)*(xgrid-xmid)
  
  f1xobs = fpm*(f1xobs)+alpha-f1min
  f1xgrid = fpm*(f1xgrid)+alpha-f1min
  
  f2xobs = fpm*f2xobs
  f2xgrid = fpm*f2xgrid
  
  return(list(fxobs=fxobs, fxgrid=fxgrid))
}


'GetRotateSf' = function(fpm,omega,psi,alpha,theta,xobs,xgrid,phixobs,phixgrid,xdelta,
						xinxgrid,xidelta,xrange,xmid,quadfacts,intsimpfacts,nbasis,nr,
						nobs,ngrid)
{	
  
  # @ Compute z2 at grid and observations @
  # @ Quadratic Form Multiplication @
  z2xobs = QuadMult(theta,phixobs,quadfacts)
  z2xgrid = QuadMult(theta,phixgrid,quadfacts)
  
  # @ Compute squish function at grid and observations @
  hxobs = SquishUp(xobs,psi,omega,nobs)
  hxgrid = SquishUp(xgrid,psi,omega,ngrid)
  
  # @ Second derivatives at xgrid and xobs @
  f2xobs = z2xobs*hxobs
  f2xgrid = z2xgrid*hxgrid
  
  # @ First derivatives at xgrid and xobs  @
  f1xgrid = InTrapCum(f2xgrid,xdelta,ngrid)
  f1xobs = IntFobs(f2xobs,f2xgrid,f1xgrid,xinxgrid,xidelta,nobs,ngrid)
  
  # @ First derivative needs to be positive. @
  f1min=min(c(0,min(f1xgrid)))  
  
  # @ Compute f at xgrid and xobs @
  fxgrid = InTrapCum(f1xgrid,xdelta,ngrid)
  fxobs = IntFobs(f1xobs,f1xgrid,fxgrid,xinxgrid,xidelta,nobs,ngrid)
  
  # @ Mean center f @
  fint = IntSimpsonfxgrid(fxgrid,xdelta,intsimpfacts,ngrid)
  fxgrid = fxgrid-fint/xrange
  fxobs = fxobs-fint/xrange	
  
  # @ Add linear term @
  fxgrid = fpm*(fxgrid)+(alpha-f1min)*(xgrid-xmid)
  fxobs = fpm*(fxobs)+(alpha-f1min)*(xobs-xmid)
  
  f1xgrid = fpm*(f1xgrid)+alpha-f1min
  f1xobs = fpm*(f1xobs)+alpha-f1min
  
  f2xgrid=fpm*f2xgrid
  f2xobs=fpm*f2xobs
  
  return(list(fxobs=fxobs, fxgrid=fxgrid))
}


'GetUf' = function(fpm,omega,psi,theta,xobs,xgrid,phixobs,phixgrid,xdelta,xinxgrid,
                 	xidelta,xrange,xmid,quadfacts,intsimpfacts,nbasis,nr,nobs,ngrid)
{	
  # @ Compute z2 at grid and observations @
  # @ Quadratic Form Multiplication @
  z2xobs = QuadMult(theta,phixobs,quadfacts)
  z2xgrid = QuadMult(theta,phixgrid,quadfacts)
  
  # @ Compute squish function at grid and observations @
  hxobs = SquishDown(xobs,psi,omega,nobs)
  hxgrid = SquishDown(xgrid,psi,omega,ngrid)
  
  # @ First derivatives at xgrid and xobs @
  f1xobs = z2xobs*hxobs
  f1xgrid = z2xgrid*hxgrid
  
  # @ Compute f at xgrid and xobs @'
  fxgrid=InTrapCum(f1xgrid,xdelta,ngrid)
  fxobs=IntFobs(f1xobs,f1xgrid,fxgrid,xinxgrid,xidelta,nobs,ngrid)
  
  # @ Mean center f @
  fint = IntSimpsonfxgrid(fxgrid,xdelta,intsimpfacts,ngrid)
  fxgrid = fxgrid-fint/xrange
  fxobs = fxobs-fint/xrange	
  
  if (fpm < 0) {
  	fxobs=-fxobs  
    fxgrid=-fxgrid
    
    f1xobs=-f1xobs
    f1xgrid=-f1xgrid
  }
  
  return(list(fxobs=fxobs, fxgrid=fxgrid))
}


#######################################
# Functions for S-shaped and U-shaped #
#######################################

#-----------------------------------------------------------------------------------------
# * SquishDown function 
# *  Decreasing logistic between -1 and 1
# *  f(x) = {1- exp[psi*(x-omega)]}/{1+exp[psi*(x-omega)]}
#-----------------------------------------------------------------------------------------
'SquishDown' = function(x,psi,omega,n)
{
  
  c = psi*(x-omega)
  
  # *******************************************
  # * Worry about underflow and INF
  # *	if c < -100, then squish = 1
  # * if c > 100,  then squish is -1
  # *******************************************
  z	= c<=-100
  c	= z*(-100-c)+c	# Change numbers <-100 to -100
  z	= c>=100
  c	= z*(100-c)+c	# Change numbers >100 to 100
  c	= exp(c)
  c	= (1-c)/(1+c)
  return(c)
 }


#-----------------------------------------------------------------------------------------
# * SquishUp function 
# *  Increasing logistic between -1 and 1
# *  f(x) = {exp[psi*(x-omega)]-1}/{exp[psi*(x-omega)]+1}
#-----------------------------------------------------------------------------------------
'SquishUp'=function(x,psi,omega,n)
{
  
  c = psi*(x-omega)
  # *******************************************
  # Worry about underflow and INF
  # *	if c < -100, then squisn = 1
  # *   if c > 100,  then squish is -1
  # *******************************************
  z	= c<=-100
  c	= z*(-100-c)+c	# Change numbers <-100 to -100
  z	= c>=100
  c	= z*(100-c)+c	# Change numbers >100 to 100
  c	= exp(c)
  c	= (c-1)/(1+c)
  return(c)
}

#-----------------------------------------------------------------------------------------
# *  INTrapCUM performs Trapizod integration 
# *  f 		= function over xgrid to be integrated 
# *  delta 	= grid size 
# *  Returns cumulative integrals over grid
# *  fint = InTrapCum(f,delta);
#-----------------------------------------------------------------------------------------

'InTrapCum' = function(f, delta, n){
  # @ Compute Trapisod approximation delta*(yi+ yi+1)/2 for each interval @
  fintmp = c(0, delta*(f[1:(n-1)]+f[2:n]) / 2)
  
  # @ Cumulative approximation @
  fint=cumsum(fintmp)

  return(fint)
}

#-----------------------------------------------------------------------------------------
# * IntFobs(hobs,fxgrid, xinxgrid,xidelta);
# *	fxgrid gives integral of h(s) at grid points xgrid
# *	We use fxgrid, which is integral of h at grid points 
# *		int_0^xi v(s) ds = int_0^x h(s)ds + int_x^xi v(s) ds
# *		= f(x) + (xi - x)*(h(xi) + h(x))/2;
# *	INPUT
# *		hobs 		= h evaluated at observations
# *		hxgrid  	= h evaluates at grid points  (finer grid than xgrid)
# *		fxgrid 		= f(x) = int_0^x h(s)ds at xgrid points
# *		xinxgrid 	= index of largest xgrid point that is less than xi
# *		xidelta		= xi - x
# *	OUTPUT
# *		f(xi)
# *
# *	fxobs = IntFobs(hobs,hxgrid, fxgrid,xinxgrid,xidelta);
#-----------------------------------------------------------------------------------------
'IntFobs' = function(hobs,hxgrid,fxgrid,xinxgrid,xidelta,nobs,ngrid)
{	
  fxobsout=fxgrid[xinxgrid]
  fxobsout=fxobsout+xidelta*(hxgrid[xinxgrid]+hobs)/2
  return(fxobsout)
}


#-----------------------------------------------------------------------------------------
# * IntSimpsonfxgrid
# *	Integratgion of fxgrid by Simpson's rule 
# *	Input
# *		k = column of fxgrid
# *	Output
# *		fint = (h/3)*(f1+4*f2+2*f3+4*f4+2*f5+...+ 2*fn-2 + 4*fn-1 + fn}
# *		n is odd
#-----------------------------------------------------------------------------------------
'IntSimpsonfxgrid' = function(fxgrid,xdelta,intsimpfacts,ngrid)
{	
  sum(fxgrid*intsimpfacts)*xdelta/3
}


################
# Cosine Basis #
################

# /*
# ****************************************************************
# * Compute Trig basis functions and their integrals on a < x < b
# *  phi_0(x) = 1/sqrt(b-a)   is constant function
# *	phi_k(x) = sqrt(2/(b-a))*cos(pi*k*(x-a)/(b-a))
# *
# *	Mean center the functions psi functions: subtract int_a^b psi_jk(x) dx/(b-a)
# ****************************************************************
# */

# @ ConstFun = phi_0(x) = ones(rows(x),1)/sqrt(b-a). phi is not mean centered. @
#fn ConstFun(x,xmin,xrange) = ones(rows(x),1)/sqrt(xrange);
'ConstFun'=function(x,xmin,xrange) 
{
  x=as.matrix(x)
  matrix(1,nr=nrow(x),nc=1)/sqrt(xrange)
}

# @ ConstFun2 = phi_0(x)^2. phi is not mean centered@
# fn ConstFun2(x,xmin,xrange) = ones(rows(x),1)/(xrange);
'ConstFun2'=function(x,xmin,xrange)
{
  x=as.matrix(x)
  matrix(1,nr=nrow(x),nc=1)/xrange
} 

# @ CosFun(x,k) = phi_k(x) = sqrt(2/(b-a))*cos(pi*k*(x-a)/(b-a)). phi is not mean centered@
# fn CosFun(x,k,xmin,xrange) = sqrt(2/(xrange))*cos(pi*k'.*(x-xmin)/(xrange));
'CosFun'=function(x,k,xmin,xrange) 
{
  tmp=sqrt(2/xrange)*cos(matrix(k,ncol=1)%*%matrix((pi*(x-xmin)/xrange),nrow=1))
  t(tmp)
}

# @ CosFun2 = phi_k(x)^2. phi is not mean centered@
# fn CosFun2(x,k,xmin,xrange) = (2/(xrange))*cos(pi*k'.*(x-xmin)/(xrange))^2;
'CosFun2'=function(x,k,xmin,xrange) 
{
  tmp=(2/(xrange))*cos(matrix(k,ncol=1)%*%matrix((pi*(x-xmin)/xrange),nrow=1))^2
  t(tmp)
}


# @ ConstCosFun(x,k) = phi_0(x)*phi_k(x). phi is not mean centered@
# fn ConstCosFun(x,k,xmin,xrange) = sqrt(2)/(xrange)*cos(pi*k'.*(x-xmin)/(xrange));
'ConstCosFun'=function(x,k,xmin,xrange)
{
  tmp=(sqrt(2)/(xrange))*cos(matrix(k,ncol=1)%*%matrix((pi*(x-xmin)/xrange),nrow=1))
  t(tmp)
}

# @ CrossProdFun(x,j,k,xmin,xrange) = phi_j(x)*phi_k(x). phi is not mean centered @
# fn CrossProdFun(x,j,k,xmin,xrange) = 
# (2/(xrange))*cos(pi*j'.*(x-xmin)/(xrange)).*cos(pi*k'.*(x-xmin)/(xrange));
'CrossProdFun'=function(x,j,k,xmin,xrange)
{
  tmpj=t(cos(matrix(j,ncol=1)%*%matrix((pi*(x-xmin)/xrange),nrow=1)))
  tmpk=t(cos(matrix(k,ncol=1)%*%matrix((pi*(x-xmin)/xrange),nrow=1)))
  (2/(xrange))*tmpj*tmpk
}

# @ IntConst2 = Int_a^x phi_0(s)^2 ds = int_a^x (1/b-a) ds then Mean Centered 	@
# fn IntConst2(x,xmin,xrange) = (x-xmin)/(xrange) -1/2;
'IntConst2'=function(x,xmin,xrange) (x-xmin)/(xrange)-1/2

# @ IntIntConst2 = int_a^x int_a^s (1/b-a) du ds then Mean Centered 	@
# fn IntIntConst2(x,xmin,xrange) = (x-xmin)^2/(2*(xrange))- (xrange)/6;
'IntIntConst2'=function(x,xmin,xrange) (x-xmin)^2/(2*(xrange))-(xrange)/6

# @ IntCos(x,k) = Int_a^x phi_0(s) phi_k(s) ds 	then Mean Centered 				@
# @ IntCos(x,k) = int_a^x sqrt(2)/(b-a)*cos(pi*j*(s-a)/(b-a))ds then Mean Centered 	@
# fn IntCos(x,k,xmin,xrange) = 
# sqrt(2)*(sin(pi*k'.*(x-xmin)/(xrange))./(pi*k')) - sqrt(2)*(1-cos(pi*k'))./((pi*k')^2);
'IntCos'=function(x,k,xmin,xrange)
{
  tmpk1=t(sin(matrix(k,ncol=1)%*%matrix((pi*(x-xmin)/xrange),nrow=1)))
  tmpk2=sqrt(2)*sweep(tmpk1,2,pi*k,FUN='/')
  tmpk3=sqrt(2)*(1-cos(pi*k))/((pi*k)^2)
  sweep(tmpk2,2,tmpk3,FUN='-')
}

# @ IntIntCos(x,k) = Int_a^x Int_a^s phi_0(t) phi_k(t) dt 	then Mean Centered 				@
# @ IntIntCos(x,k) = int_a^x int_a^s sqrt(2)/(b-a)*cos(pi*j*(u-a)/(b-a))du ds then Mean Centered @
# fn IntIntCos(x,k,xmin,xrange) = 
# sqrt(2)*(xrange)*((-cos(pi*k'.*(x-xmin)/(xrange)))./((pi*k')^2));
'IntIntCos'=function(x,k,xmin,xrange)
{
  tmpk=t(-cos(matrix(k,ncol=1)%*%matrix((pi*(x-xmin)/xrange),nrow=1)))
  sqrt(2)*(xrange)*sweep(tmpk,2,(pi*k)^2,FUN='/')
}

# @ IntCos2 = Int_a^x phi_k(s)^2 ds  			then Mean Centered 				@
# @ IntCos2 = int_a^x (2/(b-a))*cos(pi*k*(s-a)/(b-a))^2 ds	then Mean Centered @
# fn IntCos2(x,k,xmin,xrange) = sin(2*pi*(k').*(x-xmin)/(xrange))./(2*pi*(k')) + (x-xmin)/(xrange) -1/2;
'IntCos2'=function(x,k,xmin,xrange)
{
    z=(x-xmin)/xrange
    t(sin(matrix(k,ncol=1)%*%matrix(2*pi*z,nrow=1)))/(2*pi*k)+z-0.5
}
#  z=(x-xmin)/xrange)
#  tmpk1=t(sin(matrix(k,ncol=1)%*%matrix((2*pi*z,nrow=1)))
#  tmpk2=sweep(tmpk1,2,2*pi*k,FUN='/')
#  sweep(tmpk2,1,(x-xmin)/xrange-1/2,FUN='+')


# @ IntIntCos2 = Int_a^x Int_b^s phi_k(t)^2 dt ds then Mean Centered 							@
# @ IntIntCos2 = int_a^x int_a^s (2/(b-a))*(cos(pi*k*(u-a)/(b-a)))^2 du ds then Mean Centered 	@
# fn IntIntCos2(x,k,xmin,xrange) =
# (1 - cos(2*pi*(k').*(x-xmin)/(xrange))).*((xrange)./(2*pi*(k'))^2) + ((x-xmin))^2/(2*(xrange))
# - (xrange)./((2*pi*k')^2) - (xrange)/6;
'IntIntCos2'=function(x,k,xmin,xrange)
{
  tmpk1=t(cos(matrix(k,ncol=1)%*%matrix((2*pi*(x-xmin)/xrange),nrow=1)))
  tmpk2=sweep(1-tmpk1,2,xrange/((2*pi*k)^2),FUN='*')
  tmpk3=sweep(tmpk2,1,((x-xmin)^2)/(2*(xrange)),FUN='+')
  sweep(tmpk3,2,(xrange)/((2*pi*k)^2)+(xrange)/6)
}

# @ IntCosCrossProd = Int_a^x phi_k(s)*phi_j(s) ds 	then Mean Centered 								@
# @ IntCosCrossProd = int_a^x (2/(b-a)*cos(pi*j*(s-a)/(b-a))*cos(pi*k*(s-a)/(b-a)) ds then Mean Centered  @
# fn IntCosCrossProd(x,j,k,xmin,xrange) = 
# sin(pi*(j+k)'.*(x-xmin)/(xrange))./(pi*(j+k)') + sin(pi*(k-j)'.*(x-xmin)/(xrange))./(pi*(k-j)')
# - (1-cos(pi*(j+k)'))./((pi*(j+k)')^2) - (1-cos(pi*(j-k)'))./((pi*(j-k)')^2);
'IntCosCrossProd'=function(x,j,k,xmin,xrange)
{
  z=(x-xmin)/xrange
  tmpjkp1=t(sin(matrix(j+k,ncol=1)%*%matrix((pi*z),nrow=1)))
  tmpjkp2=sweep(tmpjkp1,2,pi*(j+k),FUN='/')
  tmpjkm1=t(sin(matrix(k-j,ncol=1)%*%matrix(pi*z,nrow=1)))
  tmpjkm2=sweep(tmpjkm1,2,pi*(k-j),FUN='/')
  tmpjk=tmpjkp2+tmpjkm2
  tmpstat=(1-cos(pi*(j+k)))/((pi*(j+k))^2)+(1-cos(pi*(j-k)))/((pi*(j-k))^2)
  sweep(tmpjk,2,tmpstat,FUN='-')
}

# @ IntIntCrossProd = Int_a^x Int_a^s phi_k(t)*phi_j(t) dt ds then Mean Centered @
# fn IntIntCrossProd(x,j,k,xmin,xrange) = 
# (1 - cos( pi*((j+k)').*(x-xmin)/(xrange))).*(xrange)./((pi*(j+k)')^2) +
# (1 - cos( pi*((j-k)').*(x-xmin)/(xrange))).*(xrange)./((pi*(j-k)')^2)
# - (xrange)./((pi*(j+k)')^2) - (xrange)./((pi*(j-k)')^2);
'IntIntCrossProd'=function(x,j,k,xmin,xrange)
{
  tmpjkp1=(1-t(cos(matrix(j+k,ncol=1)%*%matrix((pi*(x-xmin)/xrange),nrow=1))))*(xrange)
  tmpjkp2=sweep(tmpjkp1,2,(pi*(j+k))^2,FUN='/')
  tmpjkm1=(1-t(cos(matrix(j-k,ncol=1)%*%matrix((pi*(x-xmin)/xrange),nrow=1))))*(xrange)
  tmpjkm2=sweep(tmpjkm1,2,(pi*(j-k))^2,FUN='/')
  tmpjk=tmpjkp2+tmpjkm2
  tmpstat=xrange/((pi*(j+k))^2)+xrange/((pi*(j-k))^2)
  sweep(tmpjk,2,tmpstat,FUN='-')
}

'covariance' = function(A){
  A=as.matrix(A)
  nr=nrow(A)
  z=sapply(1:ncol(A),function(i)A[,i]-mean(A[,i]))
  return(crossprod(z)/nr)
}


'GetPhi'=function(x,xmin,xrange,phijj,phijk,phi00,phij0,nbasis){
  
  nr=(nbasis+1)*(nbasis+2)/2
  idxj0=c(0,vech(diag(1,nbasis+1))[-nr])
  idxjj=c(0,idxj0[3:nr],1)
  idx00=c(1,rep(0,nr-1))
  idxjk=rep(1,nr)-idx00-idxj0-idxjj
  idx=cbind(idx00,idxj0,idxjj,idxjk)
  idx=sapply(1:4,function(x)idx[,x]*c(1:nr))
  
  get00=phi00(x,xmin,xrange) #vector nobs
  getj0=sapply(1:nbasis,function(j)phij0(x,j,xmin,xrange))
  getjj=sapply(1:nbasis,function(j)phijj(x,j,xmin,xrange))
  getjk=NULL
  for(j in 2:nbasis){
    jk=sapply(1:(j-1),function(k)phijk(x,j,k,xmin,xrange))
    getjk=cbind(getjk,jk)
  }
  
  phi=matrix(0,nr,length(x))
  phi[idx[,1],]=get00
  phi[idx[,2],]=t(getj0)
  phi[idx[,3],]=t(getjj)
  phi[idx[,4],]=t(getjk)
  return(phi)
}

'pretrig'=function(nfree,nfunconstraint,nobs,nbasis,nint,nr,nfun,fmodel,fpm,
                   xobs,xgrid,kall,xmin,xmax,xmid,xrange){
  
  phixobsfree=array(0,dim=c(nobs,nbasis,nfree))
  phixobs=array(0,dim=c(nr,nobs,nfunconstraint))
  phixgridfree=array(0,dim=c(nint+1,nbasis,nfree))
  phixgrid=array(0,dim=c(nr,nint+1,nfunconstraint))
  
  ifree=1
  irest=1
  
  for (ifun in 1:nfun) {
    
    if (fmodel[ifun] == 1) {
      # ********************************************************************************
      # * f(x) = Phi(x)'theta
      # * vector f = Phi*theta
      # * Compute Phi matrix at observations and integration grid
      # ********************************************************************************
      phixobsfree[,,ifree] = CosFun(xobs[,ifun],kall,xmin[ifun],xrange[ifun])
      phixgridfree[,,ifree] = CosFun(xgrid[,ifun],kall,xmin[ifun],xrange[ifun])
      
      ifree=ifree+1
      
    } else {
      # ********************************************************************************
      # * Precomputes trig functions
      # *	g^c(x) = theta'Phi^c(x)*theta for c = a or b
      # *	Compute Phi^c, Vectorize it with vech, Use Vector multiplication
      # ********************************************************************************
      
      # Increasing f
      if (fmodel[ifun] == 2) {
        phixobs[,,irest]=GetPhi(xobs[,ifun],xmin[ifun],xrange[ifun],IntCos2,IntCosCrossProd,
                                IntConst2,IntCos,nbasis)
        phixgrid[,,irest]=GetPhi(xgrid[,ifun],xmin[ifun],xrange[ifun],IntCos2,IntCosCrossProd,
                                 IntConst2,IntCos,nbasis)
      }
      
      # Increasing and convex f
      if (fmodel[ifun] == 3) {
        phixobs[,,irest]=GetPhi(xobs[,ifun],xmin[ifun],xrange[ifun],IntIntCos2,IntIntCrossProd,
                                IntIntConst2,IntIntCos,nbasis)
        phixgrid[,,irest]=GetPhi(xgrid[,ifun],xmin[ifun],xrange[ifun],IntIntCos2,IntIntCrossProd, 
                                 IntIntConst2,IntIntCos,nbasis)
      }
      
      # Increasing and concave f
      if (fmodel[ifun]==4) {
        xobs2=xmin[ifun]+xmax[ifun]-xobs[,ifun]
        xgrid2=xmin[ifun]+xmax[ifun]-xgrid[,ifun]
        phixobs[,,irest]=GetPhi(xobs2,xmin[ifun],xrange[ifun],IntIntCos2,IntIntCrossProd, 
                                IntIntConst2,IntIntCos,nbasis)
        phixgrid[,,irest]=GetPhi(xgrid2,xmin[ifun],xrange[ifun],IntIntCos2,IntIntCrossProd, 
                                 IntIntConst2,IntIntCos,nbasis)
      }
      
      # "S" --> Increasing Convex to Concave f
      if (fmodel[ifun]==5) {
        phixobs[,,irest]=GetPhi(xobs[,ifun],xmin[ifun],xrange[ifun],CosFun2,CrossProdFun,
                                ConstFun2,ConstCosFun,nbasis)
        phixgrid[,,irest]=GetPhi(xgrid[,ifun],xmin[ifun],xrange[ifun],CosFun2,CrossProdFun,
                                 ConstFun2,ConstCosFun,nbasis)
      }
      
      # "Rotate S" --> Increasing Concave and Convex f
      if (fmodel[ifun]==6) {
        phixobs[,,irest]=GetPhi(xobs[,ifun],xmin[ifun],xrange[ifun],CosFun2,CrossProdFun,
                                ConstFun2,ConstCosFun,nbasis)
        phixgrid[,,irest]=GetPhi(xgrid[,ifun],xmin[ifun],xrange[ifun],CosFun2,CrossProdFun,
                                 ConstFun2,ConstCosFun,nbasis)
      }
      
      # Increasing Concave to Decreasing or inverted "U" f
      if (fmodel[ifun]==7) {
        phixobs[,,irest]=GetPhi(xobs[,ifun],xmin[ifun],xrange[ifun],CosFun2,CrossProdFun, 
                                ConstFun2,ConstCosFun,nbasis)
        phixgrid[,,irest]=GetPhi(xgrid[,ifun],xmin[ifun],xrange[ifun],CosFun2,CrossProdFun, 
                                 ConstFun2,ConstCosFun,nbasis)
      }
      irest=irest+1
    }
  }
  return(list(phixobsfree=phixobsfree,phixobs=phixobs,phixgridfree=phixgridfree,
              phixgrid=phixgrid))
}


'prefoo'=function(nobs,nbasis,nint,nr,nfun,fmodel,fpm,xobs,xgrid,xmid,xrange,
                  phixobsfree,phixobs,phixgridfree,phixgrid,
                  fxobs,fxgrid,quadfacts,theta,alpha,psi,omega,
                  xinxgrid,xidelta,intsimpfacts,xdelta)
{
  ifree=1
  irest=1
  
  for (ifun in 1:nfun) {
    
    if (fmodel[ifun] == 1) {
      foo=GetFreef(theta[2:(nbasis+1),ifun],phixobsfree[,,ifree],phixgridfree[,,ifree],
                   nbasis,nobs,nint+1)
      fxobs[,ifun]=foo$fxobs
      fxgrid[,ifun]=foo$fxgrid
      
      ifree=ifree+1
      
    } else {
      
      if (fmodel[ifun] == 2) {
        foo=GetUpf(fpm[ifun],theta[,ifun],phixobs[,,irest],phixgrid[,,irest],
                   quadfacts,nbasis,nr,nobs,nint+1)
        fxobs[,ifun]=foo$fxobs
        fxgrid[,ifun]=foo$fxgrid
      }
      
      # Increasing and convex f
      if (fmodel[ifun] == 3) {
        foo=GetConvexf(fpm[ifun],alpha[ifun],theta[,ifun],xobs[,ifun],xgrid[,ifun],
                       xmid[ifun],phixobs[,,irest],phixgrid[,,irest],
                       quadfacts,nbasis,nr,nobs,nint+1)
        fxobs[,ifun]=foo$fxobs
        fxgrid[,ifun]=foo$fxgrid
      }
      
      # Increasing and concave f
      if (fmodel[ifun]==4) {
        foo=GetConcavef(fpm[ifun],alpha[ifun],theta[,ifun],xobs[,ifun],xgrid[,ifun],
                        xmid[ifun],phixobs[,,irest],phixgrid[,,irest],
                        quadfacts,nbasis,nr,nobs,nint+1)
        fxobs[,ifun]=foo$fxobs
        fxgrid[,ifun]=foo$fxgrid
      }
      
      # "S" --> Increasing Convex to Concave f
      if (fmodel[ifun]==5) {
        foo=GetSf(fpm[ifun],omega[ifun],psi[ifun],alpha[ifun],theta[,ifun],xobs[,ifun],
                  xgrid[,ifun],phixobs[,,irest],phixgrid[,,irest],xdelta[ifun],
                  xinxgrid[,ifun],xidelta[,ifun],xrange[ifun],xmid[ifun],quadfacts,
                  intsimpfacts,nbasis,nr,nobs,nint+1)
        fxobs[,ifun]=foo$fxobs
        fxgrid[,ifun]=foo$fxgrid
        
      }
      
      # "Rotate S" --> Increasing Concave and Convex f
      if (fmodel[ifun]==6) {
        foo=GetRotateSf(fpm[ifun],omega[ifun],psi[ifun],alpha[ifun],theta[,ifun],xobs[,ifun],
                        xgrid[,ifun],phixobs[,,irest],phixgrid[,,irest],xdelta[ifun],
                        xinxgrid[,ifun],xidelta[,ifun],xrange[ifun],xmid[ifun],quadfacts,
                        intsimpfacts,nbasis,nr,nobs,nint+1)
        fxobs[,ifun]=foo$fxobs
        fxgrid[,ifun]=foo$fxgrid
      }
      
      # Increasing Concave to Decreasing or inverted "U" f
      if (fmodel[ifun]==7) {
        foo=GetUf(fpm[ifun],omega[ifun],psi[ifun],theta[,ifun],xobs[,ifun], 
                  xgrid[,ifun],phixobs[,,irest],phixgrid[,,irest],xdelta[ifun], 
                  xinxgrid[,ifun],xidelta[,ifun],xrange[ifun],xmid[ifun],quadfacts,
                  intsimpfacts,nbasis,nr,nobs,nint+1)
        fxobs[,ifun]=foo$fxobs
        fxgrid[,ifun]=foo$fxgrid
      }
      
      irest=irest+1
    }
  }
  return(list(fxobs=fxobs,fxgrid=fxgrid))
}


'Getlogg.ber'=function(fmodel,fpm,betag,betam,betas,thetag,thetam,
                          thetas,tau2g,tau2m,tau2s,gammag,gammam,
                          gammas,alphag,alpham,alphas,psig,psim,psis,
                          omegag,omegam,omegas,smcmc,nparw,nfun,nbasis,
                          iflagpsi,xmax,xmin)
{
    #@ Factor to shrink std dev @
    #if(!is.matrix(betag)) betag=as.matrix(betag)
    #if(!is.matrix(thetam)) thetam=as.matrix(thetam)
    #if(!is.matrix(thetas)) thetas=as.matrix(thetas)
    sfact=0.75
    vfact			= sfact^2	#@ Factor to shrink var    @
    
    beta_mn			= betam
    beta_sn			= sfact*betas
    beta_cov		= matrix(beta_sn,nparw,nparw)*(cor(betag))*
      matrix(beta_sn,nparw,nparw)
    beta_covi		= solve(beta_cov)
    lndetbcov		= log(det(beta_cov))
    
    theta_mn		= thetam
    theta_sn		= sfact*thetas
    theta_vn		= theta_sn^2
    theta0_mn		= thetam[1,]
    theta0_sn		= sfact*thetas[1,]
    theta0_vn		= theta0_sn^2
    # Constraint theta0 > 0
    theta0_lnpn		= numeric(nfun)
    for(i in 1:nfun){
      if (fmodel[i]>1) theta0_lnpn[i]=pnorm(-theta0_mn[i]/theta0_sn[i],
                                            lower.tail=F,log.p=T);
    }
    
    tau2_rn			= 2*(2 + (tau2m/(sfact*tau2s))^2) #@ Shrink std dev @
    tau2_sn			= tau2m*(tau2_rn - 2)
    gamma_mn		= gammam
    gamma_sn		= (sfact*gammas)
    gamma_vn		= gamma_sn^2
    gamma_lnpn		= pnorm(-gamma_mn/gamma_sn,lower.tail=F,log.p=T)
    alpha_mn 	= alpham
    alpha_sn 	= sfact*alphas
    alpha_vn 	= alpha_sn^2
    alpha_lnpn  = numeric(nfun)
    for(i in 1:nfun){
      if(fmodel[i]>2) {
        if(fpm[i]==1) { # alpha is positive
          alpha_lnpn[i]=pnorm(-alpha_mn[i]/alpha_sn[i],lower.tail=F,log.p=T);
        }else{
          alpha_lnpn[i]=pnorm(-alpha_mn[i]/alpha_sn[i],lower.tail=F,log.p=T);
        }
      }
    }
    psi_mn=psim
    psi_sn=psis*sfact
    psi_vn=psi_sn^2
    omega_mn=omegam
    omega_sn=sfact*omegas
    omega_vn=omega_sn^2
    
    psi_lnpn=numeric(nfun)
    omega_lnpn=numeric(nfun)
    for(i in 1:nfun){
      if(fmodel[i]>4.5) {
        if (iflagpsi==1){
          # Normal probability psi > 0
          if(psi_sn[i]<=0) {
            psi_sn[i]=0.001
            psi_vn[i]=psi_sn[i]^2
          }
          psi_lnpn[i]=log(1-pnorm(-psi_mn[i]/psi_sn[i]))
        }
        if(omega_sn[i]<=0){
          omega_sn[i]=0.01
          omega_vn[i]=omega_sn[i]^2
        }
        omega_lnpn[i]=log(pnorm((xmax[i]-omegam[i])/(omega_sn[i]))-pnorm((xmin[i]-omegam[i])/(omega_sn[i])))
      }
    }

    'getlogg'=function(imcmc){
      
    beta=as.matrix(betag[imcmc,])
    theta=as.matrix(thetag[,,imcmc])
    tau2=as.matrix(tau2g[imcmc,])
    gampar=as.matrix(gammag[imcmc,])
    alpha=as.matrix(alphag[imcmc,])
    psi=as.matrix(psig[imcmc,])
    omega=as.matrix(omegag[imcmc,]) 
    
    logIlik=0
    
    # @ Normal G for beta @
    residbeta=beta-beta_mn
    logIlik=logIlik-crossprod(residbeta,beta_covi%*%residbeta)/2- 
      nparw*log(2*pi)/2-lndetbcov/2
    
    for (ifun in 1:nfun){
      # @ Normal G for theta @
      residtheta=(theta[,ifun]-theta_mn[,ifun])/theta_sn[,ifun]
      if(fmodel[ifun]==1) {
        logIlik=logIlik-sum(residtheta[2:(nbasis+1)]^2)/2-
          nbasis*log(2*pi)/2-sum(log(theta_sn[2:(nbasis+1),ifun]))
      } else{ 
        logIlik=logIlik-sum(residtheta^2)/2- 
          (nbasis+1)*log(2*pi)/2-sum(log(theta_sn[,ifun]))		    
        # @ include lnP(theta0>0) @
        logIlik=logIlik-theta0_lnpn[ifun]
      }
      
      # @ IG for tau2 @
      logIlik=logIlik+LogfIG(tau2[ifun],tau2_rn[ifun],tau2_sn[ifun])
      
      # @ Truncated Normal G for gamma @
      logIlik=logIlik-((gampar[ifun]-gamma_mn[ifun])^2)/(2*gamma_vn[ifun])- 
        log(2*pi*gamma_vn[ifun])/2-gamma_lnpn[ifun]
      
      if (fmodel[ifun] > 2 && fmodel[ifun] > 7) {
        # @ Truncated normal G for alpha
        logIlik=logIlik-((alpha[ifun]-alpha_mn[ifun])^2)/(2*alpha_vn[ifun])- 
          log(2*pi*alpha_vn[ifun])-alpha_lnpn[ifun]
      }
      
      if (fmodel[ifun] >= 5) {
        if(iflagpsi == 1) {
          logIlik=logIlik-((psi[ifun]-psi_mn[ifun])^2)/(2*psi_vn[ifun])- 
            log(2*pi*psi_vn[ifun])/2-psi_lnpn[ifun]
        }
        logIlik=logIlik-((omega[ifun]-omega_mn[ifun])^2)/(2*omega_vn[ifun])- 
          log(2*pi*omega_vn[ifun])/2-omega_lnpn[ifun]
      }
    }
    return(logIlik)
  }
     logg=sapply(1:smcmc,function(imcmc)getlogg(imcmc))

  return(logg)
}


#**********************************************************************
# * Generate tau2
# *	If prior tau2 ~ IG, then generate from IG
# *	If prior tau2 ~ Exponential, 
# *		then use IG as proposal distribution, and apply Metropolis
# *	rn   = r0 + nbasis
# *	sn   = s0 + sum_k theta_k^2*exp(k*gamma)
#**********************************************************************
'gtau2'=function(nfun,theta,gamvec,iflagprior,tau2,tau2_s0,
                 tau2_u0,tau2_r0,nbasis){

  for (ifun in 1:nfun){
    
    th2gam=theta[(2:(nbasis+1)),ifun]^2/gamvec[,ifun]
    sn=sum(th2gam)
    if (iflagprior == 0) {
      # @ IG Prior is conjugate @
      tau2[ifun]=(tau2_s0+sn)/(2*rgamma(1,shape=(tau2_r0+nbasis)/2))
    } else  {
      # @ Lasso prior: use slice sampling @
      bup=tau2[ifun]-log(runif(1))/tau2_u0
      
      # @ Gamma distribution truncated below @
      tau2[ifun]=1/rndgamlo((nbasis-2)/2,sn/2,1/bup)
    }
  }
  tau2
}

#*********************************************************************
# * Generate gamma with slice sampling
#*********************************************************************
'ggamma'=function(nfun,theta,tau2,gampar,gamvec,w0,nbasis,wk,
                  zeta,kbar,kall){
  for (ifun in 1:nfun){
    
    # @ theta0 does not depend on gamma @
    ck1=theta[(2:(nbasis+1)),ifun]^2/(2*tau2[ifun])
    z=(ck1==0)
    
    # @ Worry about ck1 = 0 when thetak = 0. Doesn't appear in likelihood for gamma @
    if (sum(z) == length(z)) {
      gampar[ifun] = 1/w0
    } else {
      ckz=z*(1-ck1)+ck1 # @ ckz = 1 when theta_k = 0 @
      u1=runif(nbasis)
      bnz=gampar[ifun]+(log(ckz-log(u1)*gamvec[,ifun])-log(ckz))/kall
      if (sum(z==1)>0) {
        bnz1=numeric(sum(z!=1))
        iloop=1
        for (k in 1:nbasis){
          if (z[k] != 1) {
            bnz1[iloop]=bnz[k]
            iloop=iloop+1
          }
        }
        bmin=min(bnz1)
      } else {
        bmin=min(bnz)
      }
      u2=runif(1)
      gampar[ifun]=bmin+log(u2+(1-u2)*exp(-wk*bmin))/wk
    }
    zeta[ifun]=log(tau2[ifun])-kbar*gampar[ifun]
  }
  return(list(gampar=gampar,zeta=zeta))
}

'logitinv'=function(x)
{
  y=-log(1+exp(-x))
  exp(y)
}
