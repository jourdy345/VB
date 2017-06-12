function [v] = logh(n,p,y,X,beta,mu0,sig0,sigma2,A,B)
v = -0.5*n*log(2*pi*sigma2)-0.5/sigma2*sum((y-X*beta).^2)-0.5*logdet(sig0,'chol')-0.5*p*log(2*pi) ...
    -0.5*sum((beta-mu0).*(sig0\(beta-mu0)))-(A+1)*log(sigma2)-B/sigma2+A*log(B)-gammaln(A);
end