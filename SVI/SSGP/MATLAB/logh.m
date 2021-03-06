function [v] = logh(y,Z,alpha,gamma2,sigma2,lambda,As,Ag, ...
                    mu0l,invsig0l,n,m)
   v = -0.5*n*log(gamma2)-1/(2*gamma2)*(y-Z*alpha).'*(y-Z*alpha) ...
       -m*log(sigma2/m)-m/(2*sigma2)*(alpha.'*alpha) ...
       +log(2*As)-log(pi*(As*As+sigma2))+log(2*Ag)-log(pi*(Ag*Ag+gamma2)) ...
       +0.5*logdet(invsig0l,'chol')-0.5*(lambda-mu0l).'*invsig0l*(lambda-mu0l);
end