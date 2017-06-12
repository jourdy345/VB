% generate simulation data
n = 100;
p = 4;
beta_true = [0.3; 10; 2; 6];
X = randn(n,p);
y = X * beta_true + randn(n,1) * 1.3;

% initialize hyperparameters for prior distributions
sig0 = eye(p);
mu0  = zeros(p,1);
A = 1;
B = 1;

% variational parameters
sigq = X.'*X/n;
invsigq = invpd(sigq);
muq = zeros(p,1);
Aq = 1;
Bq = 1;

D = duplication_matrix_full(p);
D_plus = duplication_matrix_MP_inv(p);
count = 0;
lb = 0;
lbc = [];
S = 40;
nIter = 4000;

% natural parameters
gauss_lambda1 = sigq\muq;
gauss_lambda2 = -0.5*(D.'*invsigq(:));
gauss_lambda = [gauss_lambda1;gauss_lambda2];
invgam_lambda = [Aq;Bq];

muq_storage = zeros(p,nIter);
sigma2_storage = zeros(2,nIter);
d = size(gauss_lambda,1);

for i = 1:nIter
   gauss_lambda_old  = gauss_lambda;
   muq_old           = muq;
   sigq_old          = sigq;
   invgam_lambda_old = invgam_lambda;
   Aq_old            = Aq;
   Bq_old            = Bq;
   rho               = 1/(5+i);
   M_mat             = 2*D_plus*(kron(muq,eye(p)));
   S_mat             = 2*D_plus*kron(sigq,sigq)*D_plus.';
   aux               = S_mat\M_mat;
   gauss_invfisher   = [invsigq+(M_mat.'*aux) -aux.';-aux inv(S_mat)];
   invgam_fisher     = [psi(1,Aq) -1/Bq;-1/Bq Aq/Bq^2];
   
   
   % generate Monte Carlo sample
   beta_sample   = mvnrnd(muq,sigq,S);
   sigma2_sample = 1./gamrnd(Aq,1/Bq,[S,1]);
   gauss_m       = zeros(d,1);
   for j = 1:S
       beta = beta_sample(j,:).';
       gauss_m = gauss_m+[beta-muq;vech(beta*beta.'-sigq-muq*muq.')];
   end
   gauss_m       = gauss_m/S;
   invgam_m      = mean(reshape(cell2mat(arrayfun(@(sigma2) [-log(sigma2)+log(Bq)-psi(Aq);-1/sigma2+Aq/Bq],sigma2_sample,'UniformOutput',false)),2,S),2);
   gauss_grad    = zeros(d,1);
   invgam_grad   = zeros(2,1);
   lb            = 0;
   for j = 1:S
      beta        = beta_sample(j,:).';
      sigma2      = sigma2_sample(j);
      logh_v      = logh(n,p,y,X,beta,mu0,sig0,sigma2,A,B);
      logq_v      = logq(p,beta,muq,invsigq,sigma2,Aq,Bq);
      f_v         = logh_v-logq_v;
      lb          = lb+f_v;
      gauss_grad  = gauss_grad+([beta-muq;vech(beta*beta.'-sigq-muq*muq.')]-gauss_m)*f_v;
      invgam_grad = invgam_grad+([-log(sigma2)+log(Bq)-psi(Aq);-1/sigma2+Aq/Bq]-invgam_m)*f_v;
   end
   gauss_grad  = gauss_grad/(S-1);
   invgam_grad = invgam_grad/(S-1);
   gauss_lambda = (1-rho)*gauss_lambda+rho*(gauss_invfisher*gauss_grad);
   invgam_lambda = (1-rho)*invgam_lambda+rho*(invgam_fisher\invgam_grad);
   lambda1 = gauss_lambda(1:p);
   lambda2 = gauss_lambda((p+1):end);
   sigq = -0.5*(reshape(D_plus.'*lambda2,p,p))\eye(p);
   muq  = sigq*lambda1;
   if sum(eig(sigq)<0)>0
      gauss_lambda = gauss_lambda_old;
      sigq         = sigq_old;
      muq          = muq_old;
   else
       invsigq = inv(sigq);
   end
   Aq = invgam_lambda(1);
   Bq = invgam_lambda(2);
   if Aq < 0
       Aq = Aq_old;
   end
   if Bq < 0
       Bq = Bq_old;
   end
   invgam_lambda = [Aq;Bq];
   lb = lb/S;
   lbc = [lbc,lb];
   muq_storage(:,i) = muq;
   sigma2_storage(:,i) = invgam_lambda;
   disp(['LB=',num2str(lb)]);
end













