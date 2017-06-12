import numpy
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

def logh(parameters,hyperparameters,y,vphi,n,J,B,delta):
  theta = parameters[:(J+1)]
  psi   = parameters[J+1]
  sig2  = parameters[J+2]
  tau2  = parameters[J+3]
  gam   = parameters[J+4]

  # extract hyperparameters
  r_sig = hyperparameters[0]
  s_sig = hyperparameters[1]
  r_tau = hyperparameters[2]
  s_tau = hyperparameters[3]
  a     = hyperparameters[4]
  b     = hyperparameters[5]
  w0    = hyperparameters[6]
  sig0  = hyperparameters[7]

  nrows = y.size
  res   = logpdfn(theta[0],0.,np.sqrt(sig2)*sig0)+logpdfIG(sig2,r_sig/2.,s_sig/2.)+logpdfIG(tau2,r_tau/2.,s_tau/2.)+logpdfIG(psi,a,b)+np.log(w0)-w0*gam
  for i in range(nrows):
    res = res+logpdfn(y[i],delta*np.inner(theta,np.dot(vphi[i,:,:],theta)),sig2)/B*n
  for j in range(1,J+1):
    res = res+logpdfn(theta[j],0,np.sqrt(sig2)*tau2*np.exp(-j*gam))
  return res


n = 200
J = 30
x = np.random.rand(n)

# construct vphi
vphi = numpy.empty((n,J+1,J+1))
vphi[:,0,0] = x-0.5
temparg1 = numpy.outer(x,numpy.linspace(1,J,J))
temparg2 = numpy.outer(numpy.ones(n),numpy.linspace(1.,J,J))
vphi[:,0,1:] = numpy.sqrt(2.)/(numpy.pi*temparg2)*numpy.sin(numpy.pi*temparg1)-numpy.sqrt(2.)/((numpy.pi*temparg2)**2.)*(1.-numpy.cos(numpy.pi*temparg2))
vphi[:,1:,0] = vphi[:,0,1:]
temparg1 = numpy.add.outer(numpy.linspace(1.,J,J),numpy.linspace(1.,J,J))
temparg2 = numpy.multiply.outer(numpy.ones(n),temparg1)
temparg1 = numpy.multiply.outer(x,temparg1)
vphi[:,1:,1:] = numpy.sin(numpy.pi*temparg1)/(numpy.pi*temparg2)-(1.-numpy.cos(numpy.pi*temparg2))/(numpy.pi*temparg2)**2.
temparg1 = numpy.subtract.outer(numpy.linspace(1.,J,J),numpy.linspace(1.,J,J))
temparg2 = numpy.multiply.outer(numpy.ones(n),temparg1)
temparg1 = numpy.multiply.outer(x,temparg1)
vphi[:,1:,1:] = vphi[:,1:,1:]+numpy.sin(numpy.pi*temparg1)/(numpy.pi*temparg2)-(1.-numpy.cos(numpy.pi*temparg2))/(numpy.pi*temparg2)**2.
temparg1 = numpy.outer(x,numpy.linspace(1.,J,J))
temparg2 = numpy.outer(numpy.ones(n),numpy.linspace(1.,J,J))
temparg3 = numpy.outer(x,numpy.ones(J))
tempmat = numpy.sin(2.*numpy.pi*temparg1)/(2.*numpy.pi*temparg2)+temparg3-0.5
for j in range(1,J+1):
  vphi[:,j,j] = tempmat[:,(j-1)]


sigma2_m0 = 1.
sigma2_v0 = 100.
sigma2_r0 = 2.*(2.+sigma2_m0**2./sigma2_v0)
sigma2_s0 = sigma2_m0*(sigma2_r0-2.)
tau2_m0 = 1.
tau2_v0 = 100.
rtau_0 = 2.*(2.+tau2_m0**2./tau2_v0)
stau_0 = tau2_m0*(rtau_0-2.)
(r0_sig,s0_sig,r0_tau,s0_tau,a,b,w0) = (sigma2_r0,sigma2_s0,tau2_m0,tau2_v0,sigma2_r0,sigma2_s0,2.)
# m0 = numpy.zeros(p)
sig0 = 10000.
hyperparameters = [r0_sig,s0_sig,r0_tau,s0_tau,a,b,w0,sig0]
# hyperparameters = [r0_sig,s0_sig,r0_tau,s0_tau,a,b,w0,m0,sig0]
mu_post = np.zeros(J+5)
L_post  = np.eye(J+5)


grad_logh = grad(logh)
LB = 0.
LBC = np.empty(0)
S = 2
epsilon = 1.0e-16
s_k = np.ones(J+5)
rho = np.zeros(J+5)
alp = 0.3
t = 1.
# iterMax = 4000
nIter = 0
thetaContainer = np.empty(0)
cont = True
count = 0

while cont:
  z_eta = np.random.randn(S,J+5)
  grad_mu = np.zeros(J+5)
  grad_L  = np.zeros((J+5,J+5))
  for s in range(S):
    z = z_eta[s,:]
    zeta = np.dot(L_post,z)+mu_post

    theta = zeta[:(J+1)]
    psi   = np.log(np.exp(zeta[J+1])+1.)
    sig2  = np.log(np.exp(zeta[J+2])+1.)
    tau2  = np.log(np.exp(zeta[J+3])+1.)
    gam   = np.log(np.exp(zeta[J+4])+1.)
    Theta = np.hstack([theta,psi,sig2,tau2,gam])

    # subsample the dataset
    row_ind = np.random.randint(0,n,B)
    grad_Theta = grad_logh(Theta,hyperparameters,y[row_ind],vphi[row_ind,:],n,J,B,delta)
    nabla_Tinv = np.hstack([np.ones(J+1),np.exp(zeta[J+1])/(1.+np.exp(zeta[J+1])),np.exp(zeta[J+2])/(1.+np.exp(zeta[J+2])),np.exp(zeta[J+3])/(1.+np.exp(zeta[J+3])),np.exp(zeta[J+4])/(1.+np.exp(zeta[J+4]))])
    nabla_logdet = np.hstack([np.zeros(J+1),1./(1.+np.exp(zeta[J+1])),1./(1.+np.exp(zeta[J+2])),1./(1.+np.exp(zeta[J+3])),1./(1.+np.exp(zeta[J+4]))])
    grad_mu += (grad_Theta/S)*nabla_Tinv+nabla_logdet/S
    grad_L  += np.outer(grad_Theta/S*nabla_Tinv+nabla_logdet/S,z)+np.linalg.inv(L_post).T/S
  if nIter == 0:
    s_k = grad_mu*grad_mu
  else:
    s_k = alp*grad_mu*grad_mu+(1.-alp)*s_k
  rho = (nIter+1.)**(-0.5+epsilon)/(t+np.sqrt(s_k))
  mu_post += rho*grad_mu
  L_post  += np.dot(np.diag(rho),grad_L)
  LB = 0.
  Sig = np.dot(L_post,L_post.T)
  for s in range(S):
    zeta  = np.dot(L_post,np.random.randn(J+5))+mu_post
    theta = zeta[:(J+1)]
    psi   = np.log(np.exp(zeta[J+1])+1.)
    sig2  = np.log(np.exp(zeta[J+2])+1.)
    tau2  = np.log(np.exp(zeta[J+3])+1.)
    gam   = np.log(np.exp(zeta[J+4])+1.)
    Theta = np.hstack([theta,psi,sig2,tau2,gam])
    LB   += logh(Theta,hyperparameters,y,vphi,n,J,n,delta)+np.sum(zeta[(J+1):])-np.sum(np.log(1.+np.exp(zeta[(J+1):])))-logmvnpdf(zeta,mu_post,Sig)
  LB /= S
  print(nIter+1,'th LB=',LB)
  LBC = np.append(LBC,LB)
  thetaContainer = np.column_stack((thetaContainer,mu_post))
  if nIter > 200 and (nIter+1) % 100 == 0:
    diff = np.mean(LBC[-100:])-np.mean(LBC[-200:-100])
    if diff < 0. or np.abs(diff)<0.001:
      count += 1
      if count > 3:
        cont = False
    else:
      count = 0
  nIter += 1












