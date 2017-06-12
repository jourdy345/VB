function [v] = logq(p,beta,muq,invsigq,sigma2,Aq,Bq)
v = -p*0.5*log(2*pi)+0.5*logdet(invsigq,'chol')-0.5*sum((beta-muq).*(invsigq*(beta-muq)))-(Aq+1)*log(sigma2) ...
    -Bq/sigma2+Aq*log(Bq)-gammaln(Aq);
end