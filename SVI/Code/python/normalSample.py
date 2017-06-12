import autograd.numpy as np
from autograd import grad
import autograd.scipy.special as sci
from functools import partial
import matplotlib.pyplot as plt
import os
# MFVA for Normal Sample
n = 300
mu_true = 12.52
sig_true = 1.3

X = np.random.randn(n)*sig_true+mu_true

# initialize variational parameters
muq = 1.
sigq = 1.
Aq = 1.
Bq = 1.

# initialize hyperparameters
mu0 = 1.
sig0 = 1.
A0 = 1.
B0 = 1.

dif = 1.
TOL = 1.0e-10

Aq = A0 + n/2.
LB = 0.
LBold = -np.inf
LBC = np.empty(0)
while dif > TOL:
  sigq = 1/(n*Aq/Bq+1/sig0)
  muq = sigq*(n*X.mean()*Aq/Bq+mu0/sig0)
  Bq = B0+0.5*(((X-muq)*(X-muq)).sum()+n*sigq)

  LB = 0.5-0.5*n*np.log(2.*np.pi)+0.5*np.log(sigq/sig0)-((muq-mu0)*(muq-mu0)+sigq)/(2.*sig0)+A0*np.log(B0)-Aq*np.log(Bq)+gammaln(Aq)-gammaln(A0)
  dif = LB-LBold
  LBC = np.append(LBC,LB)
  print('dif:',dif)
  LBold = LB


print('Estimate for location:',muq)
print('Estimate for scale:',Bq/(Aq-1))




# Stochastic Volatility
def logpdfn(y,mu,sigma2):
  return -0.5*np.log(2.*np.pi*sigma2)-0.5/sigma2*(y-mu)*(y-mu)

def logmvnpdf(y,mu,Sig):
  p = mu.size
  (sign,logdetS) = np.linalg.slogdet(Sig)
  return -p/2.*np.log(2.*np.pi)-0.5*logdetS-0.5*np.sum((y-mu)*np.linalg.solve(Sig,y-mu))

def logh(variational_param,y,T):
  h = variational_param[:T]
  mu = variational_param[T]
  sig = variational_param[T+1]
  phi = variational_param[T+2]
  likelihood = 0.
  hprior = logpdfn(h[0],mu,sig/(1.-phi*phi))
  for i in range(T):
    likelihood = likelihood+logpdfn(y[i],0.,np.exp(h[i]))
    if i != 0:
      hprior = hprior+logpdfn(h[i],mu+phi*(h[i-1]-mu),sig)
  likelihood = likelihood+hprior-np.log(10.*np.pi)-2.*np.log(1.+mu/10.)-np.log(sig)-0.5*np.log(20.*np.pi)-np.log(sig)*np.log(sig)/200.
  return likelihood

# logh_h   = lambda h: logh(y=y,mu=mu,sig=sig,phi=phi,T=T)
# logh_mu  = lambda mu: logh(y=y,h=h,sig=sig,phi=phi,T=T)
# logh_sig = lambda sig: logh(y=y,h=h,mu=mu,phi=phi,T=T)
# logh_phi = lambda phi: logh(y=y,h=h,mu=mu,sig=sig,T=T)

# grad_h_fun   = grad(logh_h)
# grad_mu_fun  = grad(logh_mu)
# grad_sig_fun = grad(logh_sig)
# grad_phi_fun = grad(logh_phi)



# generate pseudodata
mu_true = 2.
phi_true = 0.2
sig_true = 1.3

# Uncomment to generate new pseudo-data
T = 200
h_true = np.zeros(T)
y = np.zeros(T)
h_true[0] = np.random.randn(1)[0]*np.sqrt(sig_true/(1.-phi_true*phi_true))+mu_true
y[0] = np.random.randn(1)[0]*np.exp(h_true[0]/2.)
for i in range(1,T):
  h_true[i] = np.random.randn(1)[0]*np.sqrt(sig_true)+mu_true+phi_true*(h_true[i-1]-mu_true)
  y[i] = np.random.randn(1)[0]*np.exp(h_true[i]/2.)


# os.chdir("/Users/daeyounglim/Desktop")
# y = np.loadtxt("stochasticVol.csv")
# T = 200
# initialize
h = np.random.randn(T)
sig = 2.
phi = 0.3
mu = 2.
mu_post  = np.zeros(T+3)
L_post   = np.eye(T+3)

S = 10
epsilon = 1.0e-16
tau = 1.
alpha= 0.3
s_k = np.ones(T+3)
rho = np.zeros(T+3)
LB = 0.
LBC = np.empty(0)
grad_logh = grad(logh)
iterMax = 4000
thetaContainer = np.zeros((T+3,iterMax))
for nIter in range(iterMax):

  eta = np.random.randn(S,T+3)
  grad_mu = np.zeros(T+3)
  grad_L = np.zeros((T+3,T+3))

  for s in range(S):
    z = eta[s,:]
    zeta = np.dot(L_post,z)+mu_post
    h = zeta[:T]
    mu = zeta[T]
    sig = np.exp(zeta[T+1])
    phi = np.exp(zeta[T+2])/(1.+np.exp(zeta[T+2]))
    theta = np.append(h,[mu,sig,phi])
    grad_theta = grad_logh(theta,y,T)
    # print('sig=',sig)
    # theta = np.append(np.append(zeta[:(T+1)],np.exp(zeta[T+1])),np.exp(zeta[T+2])/(1.+np.exp(zeta[T+2])))
    # grad_theta = np.append(np.append(np.append(grad_h_fun(h=theta[:T]),grad_mu_fun(mu=theta[T])),grad_sig_fun(sig=theta[T+1])),grad_phi_fun(phi=theta[T+2]))
    
    nabla_Tinv = np.append(np.ones(T+1),[theta[T+1],np.exp(zeta[T+2]/((1.+np.exp(zeta[T+2]))*(1.+np.exp(zeta[T+2]))))])
    nabla_logdet = np.append(np.zeros(T+1),[1./np.log(theta[T+1]),(1.-np.exp(zeta[T+2]))/(1.+np.exp(zeta[T+2]))])
    grad_mu = grad_mu+grad_theta*nabla_Tinv+nabla_logdet
    grad_L  = grad_L+np.outer(grad_theta*nabla_Tinv+nabla_logdet,z)+np.linalg.inv(L_post).T
  grad_mu = grad_mu/S
  grad_L = grad_L/S
  if nIter == 0:
    s_k = grad_mu*grad_mu
  else:
    s_k = alpha*grad_mu*grad_mu+(1.-alpha)*s_k
  rho = 0.1*(nIter+1.)**(-0.5+epsilon)/(tau+np.sqrt(s_k))
  mu_post = mu_post+rho*grad_mu
  L_post = L_post+np.dot(np.diag(rho),grad_L)
  LB = 0.
  Sig = np.dot(L_post,L_post.T)
  for s in range(S):
    zeta = np.dot(L_post,np.random.randn(T+3))+mu_post
    h = zeta[:T]
    mu = zeta[T]
    sig = np.exp(zeta[T+1])
    phi = np.exp(zeta[T+2])/(1.+np.exp(zeta[T+2]))
    theta = np.append(h,[mu,sig,phi])
    LB = LB+logh(theta,y,T)+np.log(sig)+np.log(np.exp(zeta[T+2])/((1.+np.exp(zeta[T+2]))*(1.+np.exp(zeta[T+2]))))-logmvnpdf(zeta,mu_post,Sig)
  LB = LB/S
  print(nIter,'th LB=',LB)
  LBC = np.append(LBC,LB)
  thetaContainer[:,nIter] = mu_post

np.savetxt("LBC.csv",LBC,delimiter=",")
np.savetxt("thetaContainer.csv",thetaContainer,delimiter=",")
np.savetxt("stochasticVol_h.csv",h_true,delimiter=",")
np.savetxt("stochasticVol_y.csv",y,delimiter=",")
np.savetxt("mu_post.csv",mu_post,delimiter=",")
np.savetxt("L_post.csv",L_post,delimiter=",")

mu_est = mu_post[T]
phi_est = np.exp(mu_post[T+2])/(1.+np.exp(mu_post[T+2]))
stochVol = mu_est+phi_est*(mu_post[:(T-1)]-mu_est)
logVol = np.log(stochVol)

plt.subplot(2,3,1)
plt.plot(LBC)
plt.ylabel("Variational LB")
plt.title("Convergence of variational lower bound",fontsize=14)
plt.ylim([-20000,-200])

# plt.subplot(2,3,2)
# plt.plot(logVol)
# plt.title("Posterior mean of log volatility h_t")

plt.subplot(2,3,2)
plt.plot(thetaContainer[T,:876])
plt.ylabel("mu")
plt.title("Convergence of mu")

plt.subplot(2,3,3)
plt.plot(thetaContainer[T+1,:876])
plt.ylabel("log(sig)")
plt.title("Convergence of log(sig)")

plt.subplot(2,3,4)
plt.plot(thetaContainer[T+2,:876])
plt.ylabel("psi")
plt.title("Convergence of psi")

plt.subplot(2,3,5)
plt.imshow(np.dot(L_post,L_post.T),cmap="hot",interpolation="nearest")

plt.show()



