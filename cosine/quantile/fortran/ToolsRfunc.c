#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

/******************************************************************************************************
 Random number generators
 ******************************************************************************************************/

void F77_SUB(rndstart)(void)
{
    GetRNGstate();
}

void F77_SUB(rndend)(void)
{
    PutRNGstate();
}

double F77_SUB(rndunif)(void)
{
    return unif_rand();
}

double F77_SUB(rndnorm)(void)
{
    return norm_rand();
}

double F77_SUB(unifrnd)(double *a, double *b)
{
    return runif(*a, *b);
}

double F77_SUB(normrnd)(double *mu, double *sigma)
{
    return rnorm(*mu, *sigma);
}

double F77_SUB(lnormrnd)(double *logmean, double *logsd)
{
    return rlnorm(*logmean, *logsd);
}

double F77_SUB(gamrnd)(double *shape, double *scale)
{
    return rgamma(*shape, *scale);
}

double F77_SUB(betarnd)(double *a, double *b)
{
    return rbeta(*a, *b);
}

double F77_SUB(exprnd)(double *scale)
{
    return rexp(*scale);
}

double F77_SUB(binomrnd)(int *n, double *p)
{
    return rbinom(*n, *p);
}

double F77_SUB(chisqrnd)(double *df)
{
    return rchisq(*df);
}

double F77_SUB(logisrnd)(double *location, double *scale)
{
    return rlogis(*location, *scale);
}

double F77_SUB(invgaussrnd)(double *mu, double *lambda)
{
    /* Random variates from inverse Gaussian distribution
     * Reference: Chhikara and Folks, The Inverse Gaussian
     * Distribution, Marcel Dekker, 1989, pp53.
     * Gordon Smyth 15 Jan 98.
     * Revised by Trevor Park 14 June 2005.
     */
    double y,r1,r2,m=*mu,la=*lambda;
    y=rchisq(1.0);
    r2=m/(2.0*la)*(2.0*la+m*y+sqrt(4.0*la*m*y+R_pow(m,2)*R_pow(y,2)));
    r1=R_pow(m,2)/r2;
    if(unif_rand() < m/(m+r1)){
        return r1;
    }else{
        return r2;
    }
}

double F77_SUB(alaplacernd)(double *location, double *scale, double *kappa)
{
    /* Return variates from asymmetric Laplace distribution 
     * Reference: VGAM R package
     */
    double x,u1,u2,m=*location,s=*scale,k=*kappa;
    u1=unif_rand();
    u2=unif_rand();
    x=m+s*log(R_pow(u1,k)/R_pow(u2,1.0/k))*M_SQRT1_2;
    return x;
}

double F77_SUB(rtgamrnd)(double *shape, double *scale, double *up)
{
    /* Random variates from right truncated gamma
     *  f(x) \propto I(x < up) x^{shape-1}exp{-x/scale}
     */
    double gup,p,x;
    gup=pgamma(*up,*shape,*scale,1,0);
    if(gup>=1.0) {
        gup=0.99999;
    }
    if(gup<=0.0) {
        gup=0.00001;
    }
    p=unif_rand()*gup;
    x=qgamma(p,*shape,*scale,1,0);
    if(x > *up) {
        return (*up - x) + x;
    }else{
        return x;
    }
}

double F77_SUB(ltgamrnd)(double *shape, double *scale, double *low)
{
    /* Random variates from left truncated gamma
     *  f(x) \propto I(x > low) x^{shape-1}exp{-x/scale}
     */
    double glow,p,x;
    glow=pgamma(*low,*shape,*scale,1,0);
    if(glow>=0.9999) {
        /* lower bound is big: use shifted exponential approximation in tail*/
        x=(*low)-log(1.0-unif_rand())*(*scale);
    }else{
        p=unif_rand()*(1.0-glow)+glow;
        x=qgamma(p,*shape,*scale,1,0);
    }
    if(x < *low) {
        return (*low - x) + x;
    }else{
        return x;
    }
}

double F77_SUB(rtnormrnd)(double *mu, double *sigma, double *up)
/* returns x ~ N(m,s^2), with x < up */
{
    double u,pcut,x;
    
    if ((*sigma)==0.0){ /* degenerate sigma=0 */
        return ((*mu)<(*up)) ? (*mu) : *(up);
    }
    pcut=pnorm(*up,*mu,*sigma,1,0);
    if (pcut < 0.0001)
        return (*up)-0.0001*(*sigma);
    u=unif_rand();
    u=u*pcut;
    x=qnorm(u,*mu,*sigma,1,0);
    
    return x;
}

double F77_SUB(ltnormrnd)(double *mu, double *sigma, double *low)
/* returns x ~ N(m,s^2), with low < x */
{
    double u,pcut,x;
    
    if ((*sigma)==0.0){ /* degenerate sigma=0 */
        return ((*mu)>(*low)) ? (*mu) : (*low);
    }
    
    pcut=pnorm(*low,*mu,*sigma,1,0);
    if (pcut>0.9999)
        return (*low)+0.0001*(*sigma);
    u=unif_rand();
    u=pcut+(1.0-pcut)*u;
    x=qnorm(u,*mu,*sigma,1,0);
    
    return x;
}

double F77_SUB(tnormrnd)(double *mu, double *sigma, double *low, double *up)
/* returns y ~ N(m,s^2), with low < y < up */
{
    double u, pleft, pright, y;
    
    pleft=pnorm(*low,*mu,*sigma,1,0);
    pright=pnorm(*up,*mu,*sigma,1,0);
    if (pleft>0.9999)
        return (*low)+0.0001*fmax2((*up)-(*low),*sigma);
    if (pright<0.0001)
        return (*up)-0.0001*fmax2((*up)-(*low),*sigma);
    u=unif_rand();
    u=pleft+(pright-pleft)*u;
    y=qnorm(u,*mu,*sigma,1,0);
    
    return y; 
}

double F77_SUB(rtlogisrnd)(double *location, double *scale, double *up)
{/* return x ~ Lo(m,s), with x < up */
    double u,pcut,x;
    if ((*scale)==0.0){ /* degenerate scale=0 */
        return ((*location)<(*up)) ? (*location) : (*up);
    }
    
    pcut=plogis(*up,*location,*scale,1,0);
    if (pcut<0.0001)
        return (*up)-0.0001*(*scale);
    u=unif_rand();
    u=u*pcut;
    x=qlogis(u,*location,*scale,1,0);
    
    return x;
}

double F77_SUB(ltlogisrnd)(double *location, double *scale, double *low)
{/* return x ~ Lo(m,s), with x > low */
    double u,pcut,x;
    if((*scale)==0.0){ /* degenerate scale=0 */
        return ((*location)>(*low)) ? (*location) : (*low);
    }
    
    pcut=plogis(*low,*location,*scale,1,0);
    if(pcut>0.9999)
        return (*low)+0.0001*(*scale);
    u=unif_rand();
    u=pcut+(1.0-pcut)*u;
    x=qlogis(u,*location,*scale,1,0);
    
    return x;
}

double F77_SUB(tlogisrnd)(double *location, double *scale, double *low, double *up)
/* returns x ~ Lo(m,s), with low < x < up */
{
    double u,pleft,pright,x;
    
    pleft=plogis(*low,*location,*scale,1,0);
    pright=plogis(*up,*location,*scale,1,0);
    if (pleft>0.9999)
        return (*low)+0.0001*fmax2((*up)-(*low),*scale);
    if (pright<0.0001)
        return (*up)-0.0001*fmax2((*up)-(*low),*scale);
    u=unif_rand();
    u=pleft+(pright-pleft)*u;
    x=qlogis(u,*location,*scale,1,0);
            
    return x;
}

double F77_SUB(snrnd)(double *xi, double *omega, double *alpha)
{
    /* This file is a component of the package 'sn' for R
     * copyright (C) 1997-2014 Adelchi Azzalini
     * returns y ~ SN(xi,omega,alpha) */
    double u1, u2, z, y;
    
    u1=norm_rand();
    u2=norm_rand();
    if (u2 > (*alpha)*u1)
        u1=-u1;
        z=u1;
        y=(*xi)+(*omega)*z;
        
        return y;
}

/******************************************************************************************************
 Cumulative density function (CDF)
 ******************************************************************************************************/

double F77_SUB(cdfcauchy)(double *x, double *location, double *scale, int *lower_tail, int *give_log)
{
	return pcauchy(*x, *location, *scale, *lower_tail, *give_log);
}

double F77_SUB(cdfnorm)(double *x, double *mu, double *sigma, int *lower_tail, int *give_log)
{
	return pnorm(*x, *mu, *sigma, *lower_tail, *give_log);
}

double F77_SUB(cdflogis)(double *x, double *location, double *scale, int *lower_tail, int *give_log)
{
	return plogis(*x, *location, *scale, *lower_tail, *give_log);
}

double F77_SUB(cdflnorm)(double *x, double *logmean, double *logsd, int *lower_tail, int *give_log)
{
	return plnorm(*x, *logmean, *logsd, *lower_tail, *give_log);
}


double F77_SUB(cdfbetas)(double *x, double *a, double *b, int *lower_tail, int *give_log)
{
	return pbeta(*x, *a, *b, *lower_tail, *give_log);
}


double F77_SUB(cdfchisq)(double *x, double *df, int *lower_tail, int *give_log)
{
	return pchisq(*x, *df, *lower_tail, *give_log);
}

double F77_SUB(cdfgamm)(double *x, double *shape, double *scale, int *lower_tail, int *give_log)
{
	return pgamma(*x, *shape, *scale, *lower_tail, *give_log);
}

double F77_SUB(cdfpoiss)(double *x, double *lambda, int *lower_tail, int *give_log)
{
	return ppois(*x, *lambda, *lower_tail, *give_log);
}

double F77_SUB(cdfalaplace)(double *x, double *location, double *scale, double *kappa)
{
    /* Return variates from asymmetric Laplace distribution
     * Reference: VGAM R package
     */
    double q=*x,m=*location,s=*scale,k=*kappa;
    double exponent,temp,p;
    if(q >= m){
        exponent=-(M_SQRT2/s)*fabs(q-m)*k;
    }else{
        exponent=-(M_SQRT2/s)*fabs(q-m)/k;
    }
    temp=exp(exponent)/(1.0+R_pow(k,2.0));
    if(q < m){
        p=R_pow(k,2.0)*temp;
    }else{
        p=1.0-temp;
    }
    return p;
}


/******************************************************************************************************
 Quantile
 ******************************************************************************************************/

double F77_SUB(invcdflogis)(double *p, double *location, double *scale, int *lower_tail, int *log_p)
{
	return qlogis(*p, *location, *scale, *lower_tail, *log_p);
}


double F77_SUB(invcdfcauchy)(double *p, double *location, double *scale, int *lower_tail, int *log_p)
{
	return qcauchy(*p, *location, *scale, *lower_tail, *log_p);
}


double F77_SUB(invcdfnorm)(double *p, double *mu, double *sigma, int *lower_tail, int *log_p)
{
	return qnorm(*p, *mu, *sigma, *lower_tail, *log_p);
}


int F77_SUB(invcdfbinom)(double *p, int *size, double *prob, int *lower_tail, int *log_p)
{
	return qbinom(*p, *size, *prob, *lower_tail, *log_p);
}


double F77_SUB(invcdfbetas)(double *p, double *a, double *b, int *lower_tail, int *log_p)
{
	return qbeta(*p, *a, *b, *lower_tail, *log_p);
}


double F77_SUB(invcdfchisq)(double *p, double *df, int *lower_tail, int *log_p)
{
	return qchisq(*p, *df, *lower_tail, *log_p);
}

double F77_SUB(invcdfpoiss)(double *p, double *lambda, int *lower_tail, int *log_p)
{
	return qpois(*p, *lambda, *lower_tail, *log_p);
}

double F77_SUB(invcdfalaplace)(double *p, double *location, double *scale, double *kappa)
{
    /* Return variates from asymmetric Laplace distribution
     * Reference: VGAM R package
     */
    double x=*p,m=*location,s=*scale,k=*kappa;
    double exponent,temp,q,inf=1.0/0.0,neginf=log(0);
    temp=R_pow(k,2.0)/(1.0+R_pow(k,2.0));
    if(x <= temp){
        exponent=x/temp;
        q=m+(s*k)*log(exponent)*M_SQRT1_2;
    }else if (x == 0.0){
        q=neginf;
    }else if (x == 1.0) {
        q=inf;
    }else{
        q=m-(s/k)*(log1p(R_pow(k,2.0))+log1p(-x))*M_SQRT1_2;
    }
    return q;
}

/******************************************************************************************************
 Probability density function (PDF)
 ******************************************************************************************************/

double F77_SUB(dbin)(double *x, double *n, double *p, int *give_log)
{
	return dbinom(*x, *n, *p, *give_log);
}


double F77_SUB(dnrm)(double *x, double *mu, double *sigma, int *give_log)
{
	return dnorm(*x, *mu, *sigma, *give_log);
}


double F77_SUB(dlnrm)(double *x, double *logmean, double *logsd, int *give_log)
{
	return dlnorm(*x, *logmean, *logsd, *give_log);
}


double F77_SUB(dlogit)(double *x, double *location, double *scale, int *give_log)
{
	return dlogis(*x, *location, *scale, *give_log);
}

double F77_SUB(dcauch)(double *x, double *location, double *scale, int *give_log)
{
	return dcauchy(*x, *location, *scale, *give_log);
}


double F77_SUB(dbet)(double *x, double *a, double *b, int *give_log)
{
	return dbeta(*x, *a, *b, *give_log);
}

double F77_SUB(dexpo)(double *x, double *scale, int *give_log)
{
    return dexp(*x, *scale, *give_log);
}

double F77_SUB(dgamm)(double *x, double *shape, double *scale, int *give_log)
{
	return dgamma(*x, *shape, *scale, *give_log);
}

double F77_SUB(dpoiss)(double *x, double *lambda, int *give_log)
{
	return dpois(*x, *lambda, *give_log);
}

double F77_SUB(dlogist)(double *x, double *location, double *scale, int *give_log)
{
    return dlogis(*x, *location, *scale, *give_log);
}

double F77_SUB(dalaplace)(double *x, double *location, double *scale, double *kappa, int *give_log)
{
    /* Return variates from asymmetric Laplace distribution
     * Reference: VGAM R package
     */
    double y=*x,m=*location,s=*scale,k=*kappa;
    double exponent,logconst,dens;
    logconst=M_LN2/2.0-log(s)+log(k)-log1p(R_pow(k,2.0));
    if(y >= m){
        exponent=-(M_SQRT2/s)*fabs(y-m)*k;
    }else{
        exponent=-(M_SQRT2/s)*fabs(y-m)/k;
    }
    dens=logconst+exponent;
    if(!*give_log)
        dens=exp(logconst+exponent);
    return dens;
}

/******************************************************************************************************
 Math functions
 ******************************************************************************************************/

double F77_SUB(gammaln)(double *x)
{
	return lgammafn(*x);
}

double F77_SUB(lbetaf)(double *a, double *b)
{
	return lbeta(*a, *b);
}

double F77_SUB(trigamm)(double *x)
{
	return trigamma(*x);
}

double F77_SUB(powerxy)(double *x, double *y)
{
    return R_pow(*x, *y);
}

double F77_SUB(powerxi)(double *x, int *i)
{
    return R_pow_di(*x, *i);
}

double F77_SUB(log1px)(double *x)
{
    return log1p(*x);
}

double F77_SUB(log1pxmx)(double *x)
{
    return log1pmx(*x);
}

/******************************************************************************************************
 Print on screen
 ******************************************************************************************************/

void F77_SUB(sprint)(int *iter, int *tot_iter, double *sec)
{
	Rprintf("MCMC draws %i of %i (CPU time: %.3f s)\n", *iter, *tot_iter, *sec);
}

void F77_SUB(progressbar)(int *c, int *t, double *sec)
{
    /* Prints a progress bar
     * INPUTS:
     * c = current iteration
     * t = total iterations
     * */
    int cit=*c,tit=*t;
    int barlength = 50;
    int nstars = barlength*cit/tit;
    Rprintf("\r|");
    for (int j=1; j<=nstars; j++){
        Rprintf("*");
    }
    Rprintf("%-*s| %*i%% (CPU time: %.3f s)", barlength - nstars, "", 3, 100*cit/tit, *sec);
}


/******************************************************************************************************
 Missing value
 ******************************************************************************************************/

int F77_SUB(isnan)(double *x)
{
    int y;
    y=R_IsNaN(*x);
    return y;
}

int F77_SUB(isna)(double *x)
{
    int y;
    y=R_IsNA(*x);
    return y;
}

int F77_SUB(ismiss)(double*x)
{
    if(R_IsNaN(*x) || R_IsNA(*x)) return 1;
    else return 0;
}
