import autograd.numpy as np
from autograd import grad
import autograd.scipy.special as sci
from functools import partial
import matplotlib.pyplot as plt
import os
import scipy.stats

T = 200
y = scipy.stats.cauchy.rvs(loc=2,scale=4,size=T)

def logpdfn(y,mu,sigma2):
  return -0.5*np.log(2.*np.pi*sigma2)-0.5/sigma2*(y-mu)*(y-mu)

def logh(theta,y,n):
  return -n*np.log(np.pi)-np.sum(np.log(1.+(y-theta)**2.))

grad_logh = grad(logh)

mu_post  = np.zeros(T+3)
L_post   = np.eye(T+3)

mu = 1.
sig = 1.

S = 10
epsilon = 1.0e-16
tau = 1.
rho = 0.
alpha = 0.3
sk = 0.
LB = 0.
LBC = np.empty(0)
iterMax = 4000
thetaContainer = np.zeros((2,iterMax))
for nIter in range(iterMax):
  eta = np.random.randn(S)
  grad_mu = 0.
  grad_sig = 0.
  for s in range(S):
    z = eta[s]
    zeta = z*sig+mu
    grad_theta = grad_logh(zeta,y,T)
    grad_mu = grad_mu+grad_theta
    grad_sig = grad_sig+grad_theta*z+1./sig
  grad_mu = grad_mu/S
  grad_sig = grad_sig/S
  if nIter == 0:
    s_k = grad_mu*grad_mu
  else:
    s_k = alpha*grad_mu*grad_mu+(1.-alpha)*s_k
  rho = 0.1*(nIter+1.)**(-0.5+epsilon)/(tau+np.sqrt(s_k))
  mu = mu+rho*grad_mu
  sig = sig+rho*grad_sig
  LB = 0.
  for s in range(S):
    zeta = np.random.randn(1)[0]*sig+mu
    LB = LB+logh(zeta,y,T)-logpdfn(zeta,mu,sig*sig)
  LB = LB/S
  print(nIter,'th LB=',LB)
  LBC = np.append(LBC,LB)
  thetaContainer[:,nIter] = np.array([mu,sig])
  print('LB:',LB)