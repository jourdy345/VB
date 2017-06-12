import theano
import numpy as np
from theano.tensor import slinalg
from scipy.special import gammaln

x = theano.tensor.dmatrix('x')
beta = theano.tensor.dvector('beta')
y = (beta*theano.dot(x,beta)).sum()
f = theano.function([beta,x],(beta*theano.dot(x,beta)).sum())
dbeta = theano.grad(y,beta)
df = theano.function([beta,x],dbeta)

dx = theano.grad(y,x)
df2 = theano.function([beta,x],dx)


lgd = theano.tensor.log(theano.tensor.nlinalg.Det()(x))
g = theano.grad(lgd,x)
gf = theano.function([x],g)
a = theano.tensor.dscalar('a')

h = theano.tensor.gammaln(a)
gh = theano.grad(h,a)
ghf = theano.function([a],gh)



sigma2 = theano.tensor.dscalar('sigma2')
y = theano.tensor.dvector('y')
X = theano.tensor.dmatrix('x')
beta = theano.tensor.dvector('beta')
mub = theano.tensor.dvector('mub')
sigb = theano.tensor.dmatrix('sigb')
A = theano.tensor.dscalar('A')
B = theano.tensor.dscalar('B')
muq = theano.tensor.dvector('muq')
sigq = theano.tensor.dmatrix('sigq')
Aq = theano.tensor.dscalar('Aq')
Bq = theano.tensor.dscalar('Bq')
n = theano.tensor.dscalar('n')

logh = -n/2.*theano.tensor.log(2.*np.pi*sigma2)-0.5/sigma2*((y-theano.dot(X,beta))*(y-theano.dot(X,beta))).sum()-0.5*((beta-mub)*slinalg.Solve()(sigb,beta-mub)).sum()-0.5*theano.tensor.log(theano.tensor.nlinalg.Det()(sigb))+A*theano.tensor.log(B)-theano.tensor.gammaln(A)-(A+1.)*theano.tensor.log(sigma2)-B/sigma2
logq = -0.5*theano.tensor.log(theano.tensor.nlinalg.Det()(sigq))-0.5*((beta-muq)*slinalg.Solve()(sigq,beta-muq)).sum()+Aq*theano.tensor.log(Bq)-theano.tensor.gammaln(Aq)-(Aq+1.)*theano.tensor.log(sigma2)-Bq/sigma2

dbeta = theano.grad(logq,beta)
dsigma2 = theano.grad(logq,sigma2)

dlogq_dbeta = theano.function([beta,muq,sigq],dbeta)
dlogq_dsigma2 = theano.function([sigma2,Aq,Bq],dsigma2)