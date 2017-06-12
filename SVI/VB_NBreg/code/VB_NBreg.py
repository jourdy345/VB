import autograd.numpy as np
from autograd import grad
import autograd.scipy.special as sci
from functools import partial
import matplotlib.pyplot as plt
import os

np.random.seed(1)
def logpdfn(y,mu,sigma2):
  return -0.5*np.log(2.*np.pi*sigma2)-0.5/sigma2*(y-mu)*(y-mu)

def logmvnpdf(y,mu,Sig):
  p = mu.size
  (sign,logdetS) = np.linalg.slogdet(Sig)
  return -p/2.*np.log(2.*np.pi)-0.5*logdetS-0.5*np.sum((y-mu)*np.linalg.solve(Sig,y-mu))

def logpdfIG(x,A,B):
  return A*np.log(B)-sci.gammaln(A)-(A+1.)*np.log(x)-B/x

# n: length of y
# p: length of beta
# J: length of theta (truncation point of Gaussian process)
def logh(variational_param,hyperprior,y,vphi,W,p,J,n):
  # extract parameters
  beta  = variational_param[:p]
  theta = variational_param[p:(p+J)]
  psi   = variational_param[p+J]
  sig2  = variational_param[p+J+1]
  tau2  = variational_param[p+J+2]
  gam   = variational_param[p+J+3]

  # extract hyperparameters
  r_sig = hyperprior[0]
  s_sig = hyperprior[1]
  r_tau = hyperprior[2]
  s_tau = hyperprior[3]
  a     = hyperprior[4]
  b     = hyperprior[5]
  w0    = hyperprior[6]
  m0    = hyperprior[7]
  sig0  = hyperprior[8]

  res   = n*sci.gammaln(psi)+n*psi*np.log(psi)+np.inner(y,np.dot(W,beta))+np.inner(y,np.dot(vphi,theta))+logpdfIG(sig2,r_sig/2.,s_sig/2.)+logpdfIG(tau2,r_tau/2.,s_tau/2.)+logpdfIG(psi,a,b)+np.log(w0)-w0*gam+logmvnpdf(beta,m0,sig2*sig0)
  for i in range(n):
    # res += sci.gammaln(y[i]+psi)-sci.gammaln(y[i]+1.)-(y[i]+psi)*np.log(psi+np.exp(np.inner(W[i,:],beta)+np.inner(vphi[i,:],theta)))
    res = res+sci.gammaln(y[i]+psi)-(y[i]+psi)*np.log(psi+np.exp(np.inner(W[i,:],beta)+np.inner(vphi[i,:],theta)))
  for j in range(J):
    res = res+logpdfn(theta[j],0.,sig2*tau2*np.exp(-(j+1.)*gam))
  return res

# def f(x):
#   return np.exp(6.*x-3.)
def f(x):
  return 16.*x**2.-4.*np.cos(2.*np.pi*x)/(np.pi**2.)-np.cos(4.*np.pi*x)/(np.pi**2.)-32.*np.cos(3.*np.pi*x)/(9.*np.pi**2.)-32.*np.cos(np.pi*x)/(np.pi**2.)+365./(9.*np.pi**2.)

n = 200
p = 2
x = np.random.rand(n)
fx = f(x)
xmin = np.min(x)
xmax = np.max(x)
x = (x-xmin)/(xmax-xmin)
y = np.empty(n)

J = 10
psi = 2.4
beta_true = np.array([1.4,-2.5])
W = np.random.rand(n,p)
W_beta = np.dot(W,beta_true)
vphi = np.sqrt(2.)*np.cos(np.outer(x,np.pi*np.array(range(1,(J+1)))))
for k in range(n):
  tmp = np.random.gamma(shape=psi,scale=np.exp(fx[k]+W_beta[k])/psi,size=1)[0]
  y[k] = np.random.poisson(lam=tmp,size=1)[0]


(r0_sig,s0_sig,r0_tau,s0_tau,a,b,w0) = (10.,10.,10.,10.,10.,10.,2.)
m0 = np.zeros(p)
sig0 = np.eye(p)
hyperprior = np.array([r0_sig,s0_sig,r0_tau,s0_tau,a,b,w0,m0,sig0])
mu_post = np.zeros(p+J+4)
L_post  = np.eye(p+J+4)


grad_logh = grad(logh)
LB = 0.
LBC = np.empty(0)
S = 20
epsilon = 1.0e-16
s_k = np.ones(p+J+4)
rho = np.zeros(p+J+4)
alp = 0.3
t = 1.
iterMax = 1000
thetaContainer = np.zeros((p+J+4,iterMax))
for nIter in range(iterMax):
  z_eta = np.random.randn(S,p+J+4)
  grad_mu = np.zeros(p+J+4)
  grad_L = np.zeros((p+J+4,p+J+4))
  for s in range(S):
    z = z_eta[s,:]
    zeta = np.dot(L_post,z)+mu_post
    beta  = zeta[:p]
    theta = zeta[p:(p+J)]
    psi   = np.log(np.exp(zeta[p+J])+1.)
    sig2  = np.log(np.exp(zeta[p+J+1])+1.)
    tau2  = np.log(np.exp(zeta[p+J+2])+1.)
    gam   = np.log(np.exp(zeta[p+J+3])+1.)
    Theta = np.hstack([beta,theta,psi,sig2,tau2,gam])
    grad_Theta = grad_logh(Theta,hyperprior,y,vphi,W,p,J,n)
    nabla_Tinv = np.hstack([np.ones(p+J),np.exp(zeta[p+J])/(1.+np.exp(zeta[p+J])),np.exp(zeta[p+J+1])/(1.+np.exp(zeta[p+J+1])),np.exp(zeta[p+J+2])/(1.+np.exp(zeta[p+J+2])),np.exp(zeta[p+J+3])/(1.+np.exp(zeta[p+J+3]))])
    nabla_logdet = np.hstack([np.zeros(p+J),1./(1.+np.exp(zeta[p+J])),1./(1.+np.exp(zeta[p+J+1])),1./(1.+np.exp(zeta[p+J+2])),1./(1.+np.exp(zeta[p+J+3]))])
    grad_mu += (grad_Theta/S)*nabla_Tinv+nabla_logdet/S
    grad_L  += np.outer(grad_Theta/S*nabla_Tinv+nabla_logdet/S,z)+np.linalg.inv(L_post).T/S
  if nIter == 0:
    s_k = grad_mu*grad_mu
  else:
    s_k = alp*grad_mu*grad_mu+(1.-alp)*s_k
  rho = 0.1*(nIter+1.)**(-0.5+epsilon)/(t+np.sqrt(s_k))
  mu_post += rho*grad_mu
  L_post  += np.dot(np.diag(rho),grad_L)
  LB = 0.
  Sig = np.dot(L_post,L_post.T)
  for s in range(S):
    zeta  = np.dot(L_post,np.random.randn(p+J+4))+mu_post
    beta  = zeta[:p]
    theta = zeta[p:(p+J)]
    psi   = np.log(np.exp(zeta[p+J])+1.)
    sig2  = np.log(np.exp(zeta[p+J+1])+1.)
    tau2  = np.log(np.exp(zeta[p+J+2])+1.)
    gam   = np.log(np.exp(zeta[p+J+3])+1.)
    Theta = np.hstack([beta,theta,psi,sig2,tau2,gam])
    LB   += logh(Theta,hyperprior,y,vphi,W,p,J,n)+np.sum(zeta[(p+J):])-np.sum(np.log(1.+np.exp(zeta[(p+J):])))-logmvnpdf(zeta,mu_post,Sig)
  LB /= S
  print(nIter,'th LB=',LB)
  LBC = np.append(LBC,LB)
  thetaContainer[:,nIter] = mu_post

os.chdir("/Users/daeyounglim/Desktop/Github/VB/SVI/VB_NBreg/code/")
np.savetxt("LBC.csv",LBC,delimiter=",")
np.savetxt("thetaContainer.csv",thetaContainer,delimiter=",")


LBC = np.loadtxt("LBC.csv",delimiter=",")
thetaContainer = np.loadtxt("thetaContainer.csv",delimiter=",")


beta  = thetaContainer[:p,-1]
for i in range(1,17):
  plt.subplot(4,4,i)
  plt.plot(thetaContainer[(i-1),:])
plt.show()

