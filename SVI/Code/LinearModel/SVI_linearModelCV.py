import numpy as np
from scipy.special import gammaln,polygamma,digamma

def duplication_matrix(n):
  k = np.zeros((n,n))
  index = np.linspace(1,n*(n+1)/2,n*(n+1)/2)
  for j in range(n):
    for i in range(j,n):
      if j == 0:
        k[i,j] = index[i]
      else:
        k[i,j] = index[n*j-(j-1)*j//2+i-j]
  k[np.triu_indices(n,1)] = k[np.tril_indices(n,-1)]
  tmpVec = np.ravel(k)
  t = tmpVec.size
  s = index.size
  res = np.zeros((t,s))
  for k in range(t):
    for u in range(s):
      if tmpVec[k] == index[u]:
        res[k,u] = 1
      else:
        res[k,u] = 0
  return res

def vech(X):
  (n,_) = X.shape
  res = np.zeros(n*(n+1)//2)
  for j in range(n):
    for i in range(j,n):
      res[n*j-(j-1)*j//2+i-j] = X[i,j]
  return res


def logh(n,p,y,X,beta,mu0,sig0,sigma2,A,B):
  y_Xb = y-X.dot(beta)
  (sign,logdet) = np.linalg.slogdet(sig0)
  v = -0.5*n*np.log(2.*sigma2*np.pi)-0.5/sigma2*np.inner(y_Xb,y_Xb)-0.5*logdet-0.5*p*np.log(2*np.pi)-0.5*np.inner(beta-mu0,np.linalg.solve(sig0,beta-mu0))-(A+1)*np.log(sigma2)-B/sigma2+A*np.log(B)-gammaln(A)
  return v

def logq(p,beta,muq,invsigq,sigma2,Aq,Bq):
  (sign,logdet) = np.linalg.slogdet(invsigq)
  invsigq_beta_muq = invsigq.dot(beta-muq)
  v = -0.5*p*np.log(2.*np.pi)+0.5*logdet-0.5*np.inner(beta-muq,invsigq_beta_muq)-(Aq+1)*np.log(sigma2)-Bq/sigma2+Aq*np.log(Bq)-gammaln(Aq)
  return v


n = 100
p = 4
beta_true = np.array([0.3,10.,2.,6.])
X = np.random.randn(100,4)
y = X.dot(beta_true)+np.random.randn(100)*1.3

# initialize hyperparameters for prior distributions
sig0 = np.eye(p)
mu0  = np.zeros(p)
(A,B,Aq,Bq) = (1.,1.,1.,1.)

D      = duplication_matrix(p)
D_plus = np.linalg.pinv(D)
count  = 0
lb     = 0.
lbc    = np.empty((0,0))
S      = 40
nIter  = 10000

# variational parameters
sigq = np.random.randn(p,p)
sigq = sigq.T.dot(sigq)
invsigq = np.linalg.inv(sigq)
muq = np.zeros(p)


# natural parameters
gauss_lambda1 = np.linalg.solve(sigq,muq)
gauss_lambda2 = -0.5*D.T.dot(np.ravel(invsigq))
gauss_lambda  = np.append(gauss_lambda1,gauss_lambda2)
invgam_lambda = np.array([Aq,Bq])

muq_storage    = np.empty((p,nIter))
sigma2_storage = np.empty((2,nIter))

# old parameters
gauss_lambda_old  = np.empty_like(gauss_lambda)
muq_old           = np.empty_like(muq)
sigq_old          = np.empty_like(sigq)
invgam_lambda_old = np.empty_like(invgam_lambda)



for i in range(nIter):
  np.copyto(gauss_lambda_old,gauss_lambda)
  np.copyto(muq_old,muq)
  np.copyto(sigq_old,sigq)
  np.copyto(invgam_lambda_old,invgam_lambda)
  (Aq_old,Bq_old) = (Aq,Bq)
  rho = 1./(6.+i)
  M_mat = 2.*D_plus.dot(np.kron(muq,np.eye(p)).T)
  S_mat = 2.*D_plus.dot(np.kron(sigq,sigq).dot(D_plus.T))
  aux   = np.linalg.solve(S_mat,M_mat)
  gauss_invfisher = np.vstack((np.hstack((invsigq+M_mat.T.dot(aux),-aux.T)),np.hstack((-aux,np.linalg.inv(S_mat)))))
  invgam_fisher   = np.array([polygamma(1,Aq),-1./Bq,-1./Bq,Aq/(Bq*Bq)]).reshape(2,2)


  beta_sample = np.random.multivariate_normal(muq,sigq,S)
  sigma2_sample = 1./np.random.gamma(shape=Aq,scale=1./Bq,size=S)
  gauss_m       = np.array(list(map(lambda beta: np.append(beta-muq,vech(np.outer(beta,beta)-sigq-np.outer(muq,muq))),beta_sample))).mean(axis=0)
  invgam_m      = np.array(list(map(lambda sigma2: np.array([-np.log(sigma2)+np.log(Bq)-digamma(Aq),-1./sigma2+Aq/Bq]),sigma2_sample))).mean(axis=0)
  gauss_grad    = np.empty_like(gauss_lambda)
  invgam_grad   = np.empty_like(invgam_lambda)
  lb            = 0.
  for j in range(S):
    beta = beta_sample[j,:]
    sigma2 = sigma2_sample[j]
    logh_v = logh(n,p,y,X,beta,mu0,sig0,sigma2,A,B)
    logq_v = logq(p,beta,muq,invsigq,sigma2,Aq,Bq)
    f_v    = logh_v-logq_v
    lb    += f_v
    gauss_grad  += (np.append(beta-muq,vech(np.outer(beta,beta)-sigq-np.outer(muq,muq)))-gauss_m)*f_v
    invgam_grad += (np.array([-np.log(sigma2)+np.log(Bq)-digamma(Aq),-1./sigma2+Aq/Bq])-invgam_m)*f_v
  gauss_grad   /= (S-1.)
  invgam_grad  /= (S-1.)
  gauss_lambda  = (1.-rho)*gauss_lambda+rho*gauss_invfisher.dot(gauss_grad)
  invgam_lambda = (1.-rho)*invgam_lambda+rho*np.linalg.solve(invgam_fisher,invgam_grad)
  gauss_lambda1 = gauss_lambda[:p]
  gauss_lambda2 = gauss_lambda[p:]
  sigq          = -0.5*np.linalg.inv(D_plus.T.dot(gauss_lambda2).reshape(p,p))
  muq           = sigq.dot(gauss_lambda1)
  try:
    pd = np.linalg.cholesky(sigq)
    invsigq = np.linalg.inv(sigq)
  except np.linalg.LinAlgError:
    np.copyto(gauss_lambda,gauss_lambda_old)
    np.copyto(muq,muq_old)
    np.copyto(sigq,sigq_old)
  Aq = invgam_lambda[0]
  Bq = invgam_lambda[1]
  if Aq<0.:
    Aq = Aq_old
  if Bq<0.:
    Bq = Bq_old
  invgam_lambda = np.array([Aq,Bq])
  lbc = np.append(lbc,lb)
  muq_storage[:,i] = muq
  sigma2_storage[:,i] = invgam_lambda
  print("LB=",lb)










