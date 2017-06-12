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

# def logh(parameters,hyperparameters,y,vphi,W,n,p,J:
#   # extract parameters
#   beta  = parameters[:p]
#   theta = parameters[p:(J)]
#   psi   = parameters[J]
#   sig2  = parameters[J+1]
#   tau2  = parameters[J+2]
#   gam   = parameters[J+3]

#   # extract hyperparameters
#   r_sig = hyperparameters[0]
#   s_sig = hyperparameters[1]
#   r_tau = hyperparameters[2]
#   s_tau = hyperparameters[3]
#   a     = hyperparameters[4]
#   b     = hyperparameters[5]
#   w0    = hyperparameters[6]
#   m0    = hyperparameters[7]
#   sig0  = hyperparameters[8]


#   res = np.inner(y,np.dot(W,beta))+np.inner(y,np.dot(vphi,theta))-np.sum(np.log(1.+np.exp(np.dot(W,beta)+np.dot(vphi,theta))))+logmvnpdf(beta,m0,sig2*sig0)+logpdfIG(sig2,r_sig/2.,s_sig/2.)+logpdfIG(tau2,r_tau/2.,s_tau/2.)+logpdfIG(psi,a,b)+np.log(w0)-w0*gam
#   for j in range(J):
#     res = res+logpdfn(theta[j],0,sig2*tau2*np.exp(-(j+1.)*gam))
#   return res

def logh(parameters,hyperparameters,y,vphi,n,J):
  # extract parameters
  theta = parameters[:J]
  psi   = parameters[J]
  sig2  = parameters[J+1]
  tau2  = parameters[J+2]
  gam   = parameters[J+3]

  # extract hyperparameters
  r_sig = hyperparameters[0]
  s_sig = hyperparameters[1]
  r_tau = hyperparameters[2]
  s_tau = hyperparameters[3]
  a     = hyperparameters[4]
  b     = hyperparameters[5]
  w0    = hyperparameters[6]


  res = logmvnpdf(y,np.dot(vphi,theta),sig2*np.eye(n))+logpdfIG(sig2,r_sig/2.,s_sig/2.)+logpdfIG(tau2,r_tau/2.,s_tau/2.)+logpdfIG(psi,a,b)+np.log(w0)-w0*gam
  for j in range(J):
    res = res+logpdfn(theta[j],0,sig2*tau2*np.exp(-(j+1.)*gam))
  return res


# def f(x):
#   return 16.*x**2.-4.*np.cos(2.*np.pi*x)/(np.pi**2.)-np.cos(4.*np.pi*x)/(np.pi**2.)-32.*np.cos(3.*np.pi*x)/(9.*np.pi**2.)-32.*np.cos(np.pi*x)/(np.pi**2.)+365./(9.*np.pi**2.)-10.
def f(x):
  return 2.*np.pi*x+np.sin(2.*np.pi*x)

n = 200
# p = 4
x = np.random.rand(n)
fx = f(x)
sig2_true = 3.4
xmin = np.min(x)
xmax = np.max(x)
x = (x-xmin)/(xmax-xmin)
y = fx+np.random.randn(n)*np.sqrt(sig2_true)


J = 20
psi = 2.4
# beta_true = np.array([1.4,-2.6,0.3,9.25])
# W = np.random.rand(n,p)
# W_beta = np.dot(W,beta_true)
vphi = np.sqrt(2.)*np.cos(np.outer(x,np.pi*np.array(range(1,(J+1)))))
# mean_parameter = fx+W_beta
# for k in range(n):
#   y[k] = 1 if fx[k]+np.random.logistic(size=1)[0] > 0. else 0


(r0_sig,s0_sig,r0_tau,s0_tau,a,b,w0) = (10.,10.,10.,10.,10.,10.,2.)
# m0 = np.zeros(p)
# sig0 = np.eye(p)
hyperparameters = [r0_sig,s0_sig,r0_tau,s0_tau,a,b,w0]
# hyperparameters = [r0_sig,s0_sig,r0_tau,s0_tau,a,b,w0,m0,sig0]
mu_post = np.zeros(J+4)
L_post  = np.eye(J+4)


grad_logh = grad(logh)
LB = 0.
LBC = np.empty(0)
S = 20
epsilon = 1.0e-16
s_k = np.ones(J+4)
rho = np.zeros(J+4)
alp = 0.3
t = 1.
iterMax = 4000
thetaContainer = np.zeros((J+4,iterMax))
for nIter in range(iterMax):
  z_eta = np.random.randn(S,J+4)
  grad_mu = np.zeros(J+4)
  grad_L = np.zeros((J+4,J+4))
  for s in range(S):
    z = z_eta[s,:]
    zeta = np.dot(L_post,z)+mu_post
    # beta  = zeta[:p]
    # theta = zeta[p:(J)]
    # psi   = np.log(np.exp(zeta[J])+1.)
    # sig2  = np.log(np.exp(zeta[J+1])+1.)
    # tau2  = np.log(np.exp(zeta[J+2])+1.)
    # gam   = np.log(np.exp(zeta[J+3])+1.)

    theta = zeta[:J]
    psi   = np.log(np.exp(zeta[J])+1.)
    sig2  = np.log(np.exp(zeta[J+1])+1.)
    tau2  = np.log(np.exp(zeta[J+2])+1.)
    gam   = np.log(np.exp(zeta[J+3])+1.)
    Theta = np.hstack([theta,psi,sig2,tau2,gam])
    grad_Theta = grad_logh(Theta,hyperparameters,y,vphi,n,J)
    nabla_Tinv = np.hstack([np.ones(J),np.exp(zeta[J])/(1.+np.exp(zeta[J])),np.exp(zeta[J+1])/(1.+np.exp(zeta[J+1])),np.exp(zeta[J+2])/(1.+np.exp(zeta[J+2])),np.exp(zeta[J+3])/(1.+np.exp(zeta[J+3]))])
    nabla_logdet = np.hstack([np.zeros(J),1./(1.+np.exp(zeta[J])),1./(1.+np.exp(zeta[J+1])),1./(1.+np.exp(zeta[J+2])),1./(1.+np.exp(zeta[J+3]))])
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
    zeta  = np.dot(L_post,np.random.randn(J+4))+mu_post
    theta = zeta[:J]
    psi   = np.log(np.exp(zeta[J])+1.)
    sig2  = np.log(np.exp(zeta[J+1])+1.)
    tau2  = np.log(np.exp(zeta[J+2])+1.)
    gam   = np.log(np.exp(zeta[J+3])+1.)
    Theta = np.hstack([theta,psi,sig2,tau2,gam])
    LB   += logh(Theta,hyperparameters,y,vphi,n,J)+np.sum(zeta[(J):])-np.sum(np.log(1.+np.exp(zeta[J:])))-logmvnpdf(zeta,mu_post,Sig)
  LB /= S
  print(nIter+1,'th LB=',LB)
  LBC = np.append(LBC,LB)
  thetaContainer[:,nIter] = mu_post


os.chdir("/Users/daeyounglim/Desktop/Github/VB/SVI/VB_logit/code/")
np.savetxt("LBC.csv",LBC,delimiter=",")
np.savetxt("thetaContainer.csv",thetaContainer,delimiter=",")

for i in range(24):
  plt.subplot(5,5,i+1)
  plt.plot(thetaContainer[i,:])
theta_est = mu_post[:J]
# f_est = np.dot(W,beta_est)+np.dot(vphi,theta_est)
f_est = np.dot(vphi,theta_est)
plt.subplot(5,5,25)
plt.plot(np.sort(x),f_est[np.argsort(x)])
plt.plot(np.sort(x),fx[np.argsort(x)],'r--')
plt.show()

# beta_est = mu_post[:p]
# theta_est = mu_post[p:J]
