% n = 100;
% d = 8;
m = 10;

load ../R/PendulumData/X_tr.txt;
load ../R/PendulumData/T_tr.txt;
y = T_tr;
[n,d] = size(X_tr);
% sig2_true = 1.3;
% x = randn(n,d);
% beta_true = randn(d,1);
% y = cos(x)*beta_true;

S_sample = mvnrnd(zeros(d,1),eye(d),m);


nIter = 7000;
S = 40;
lbc = [];
lb  = 0;
D_alpha = duplication_matrix_full(2*m);
D_plus_alpha = duplication_matrix_MP_inv(2*m);
D_lambda = duplication_matrix_full(d);
D_plus_lambda = duplication_matrix_MP_inv(d);

% hyperparameters
As = 20;
Ag = 20;
mu0l = zeros(d,1);
sig0l = eye(d);
invsig0l = invpd(sig0l);

% variational parameters
muaq  = zeros(2*m,1);
sigaq = eye(2*m);
invsigaq = invpd(sigaq);
mulq  = zeros(d,1);
siglq = eye(d);
invsiglq = invpd(siglq);
as = 1;
bs = 1;
ag = 1;
bg = 1;
aux_alpha  = eye(2*m)/sigaq;
aux_lambda = eye(d)/siglq;

alpha_lambda1 = sigaq\muaq;
alpha_lambda2 = -0.5*D_alpha'*aux_alpha(:);
lambda_lambda1 = siglq\mulq;
lambda_lambda2 = -0.5*D_lambda'*aux_lambda(:);
alpha_lambda = [alpha_lambda1;alpha_lambda2];
lambda_lambda = [lambda_lambda1;lambda_lambda2];
sigma_lambda = [as;bs];
gamma_lambda = [ag;bg];

for i = 1:nIter
    alpha_lambda_old = alpha_lambda;
    lambda_lambda_old = lambda_lambda;
    sigma_lambda_old = sigma_lambda;
    gamma_lambda_old = gamma_lambda;
    muaq_old         = muaq;
    sigaq_old        = sigaq;
    mulq_old         = mulq;
    siglq_old        = siglq;
    as_old           = as;
    bs_old           = bs;
    ag_old           = ag;
    bg_old           = bg;
    rho = 1/(5+i);
    M_alpha = 2*D_plus_alpha*kron(muaq,eye(2*m));
    S_alpha = 2*D_plus_alpha*kron(sigaq,sigaq)*D_plus_alpha';
    aux_alpha = S_alpha\M_alpha;
    alpha_invfisher = [eye(2*m)/sigaq+M_alpha'*aux_alpha -aux_alpha';-aux_alpha eye(size(S_alpha))/S_alpha];
    M_lambda = 2*D_plus_lambda*kron(mulq,eye(d));
    S_lambda = 2*D_plus_lambda*kron(siglq,siglq)*D_plus_lambda';
    aux_lambda = S_lambda\M_lambda;
    lambda_invfisher = [eye(d)/siglq+M_lambda'*aux_lambda -aux_lambda';-aux_lambda eye(size(S_lambda))/S_lambda];
    sigma_fisher = [psi(1,as) -1/bs;-1/bs as/(bs*bs)];
    gamma_fisher = [psi(1,ag) -1/bg;-1/bg ag/(bg*bg)];
    
    alpha_sample = mvnrnd(muaq,sigaq,S);
    lambda_sample = mvnrnd(mulq,siglq,S);
    sigma_sample = 1./gamrnd(as,1/bs,[S,1]);
    gamma_sample = 1./gamrnd(ag,1/bg,[S,1]);
    
    % control variate
    alpha_m = zeros(S,2*m+2*m*(2*m+1)/2);
    lambda_m = zeros(S,d+d*(d+1)/2);
    sigma_m = zeros(S,2);
    gamma_m = zeros(S,2);
    for k = 1:S
       alpha = alpha_sample(k,:)';
       lambda = lambda_sample(k,:)';
       sigma2 = sigma_sample(k);
       gamma2 = gamma_sample(k);
       alpha_m(k,:) = [alpha-muaq;vech(alpha*alpha'-sigaq-muaq*muaq')]';
       lambda_m(k,:) = [lambda-mulq;vech(lambda*lambda'-siglq-mulq*mulq')]';
       sigma_m(k,:) = [-log(sigma2)+log(bs)-psi(as);-1/sigma2+as/bs]';
       gamma_m(k,:) = [-log(gamma2)+log(bg)-psi(ag);-1/gamma2+ag/bg]';
    end
    alpha_m  = mean(alpha_m)';
    lambda_m = mean(lambda_m)';
    sigma_m  = mean(sigma_m)';
    gamma_m  = mean(gamma_m)';
   
    % gradient & lower bound
    alpha_grad = zeros(2*m+2*m*(2*m+1)/2,1);
    lambda_grad = zeros(d+d*(d+1)/2,1);
    sigma_grad = zeros(2,1);
    gamma_grad = zeros(2,1);
    lb         = 0;
    f_v_c      = zeros(1,S);
    for j = 1:S
       alpha = alpha_sample(j,:)';
       lambda = lambda_sample(j,:)';
       sigma2 = sigma_sample(j);
       gamma2 = gamma_sample(j);
       Z      = Zmat(X_tr,S_sample,lambda,n,m,d);
       logh_v = logh(y,Z,alpha,gamma2,sigma2,lambda,As,Ag,mu0l,invsig0l,n,m);
       logq_v = logq(alpha,lambda,sigma2,gamma2,muaq,invsigaq,mulq,invsiglq,as,bs,ag,bg,m,d);
       f_v    = logh_v-logq_v;
       f_v_c(j)  = logq_v;
       lb     = lb+f_v;
       alpha_grad = alpha_grad+([alpha-muaq;vech(alpha*alpha'-sigaq-muaq*muaq')]-alpha_m)*f_v;
       lambda_grad = lambda_grad+([lambda-mulq;vech(lambda*lambda'-siglq-mulq*mulq')]-lambda_m)*f_v;
       sigma_grad = sigma_grad+([-log(sigma2)+log(bs)-psi(as);-1/sigma2+as/bs]-sigma_m)*f_v;
       gamma_grad = gamma_grad+([-log(gamma2)+log(bg)-psi(ag);-1/gamma2+ag/bg]-gamma_m)*f_v;
    end
    alpha_grad = alpha_grad/(S-1);
    lambda_grad = lambda_grad/(S-1);
    sigma_grad = sigma_grad/(S-1);
    gamma_grad = gamma_grad/(S-1);
    alpha_lambda = (1-rho)*alpha_lambda+rho*alpha_invfisher*alpha_grad;
    lambda_lambda = (1-rho)*lambda_lambda+rho*lambda_invfisher*lambda_grad;
    sigma_lambda = (1-rho)*sigma_lambda+rho*sigma_fisher\sigma_grad;
    gamma_lambda = (1-rho)*gamma_lambda+rho*gamma_fisher\gamma_grad;
    
    % natural parameters -> ordinary parameters
    alpha_lambda1 = alpha_lambda(1:(2*m));
    alpha_lambda2 = alpha_lambda((2*m+1):end);
    sigaq         = -0.5*eye(2*m)/reshape(D_plus_alpha'*alpha_lambda2,2*m,2*m);
    muaq          = sigaq*alpha_lambda1;
    
    lambda_lambda1 = lambda_lambda(1:d);
    lambda_lambda2 = lambda_lambda((d+1):end);
    siglq          = -0.5*eye(d)/reshape(D_plus_lambda'*lambda_lambda2,d,d);
    mulq           = siglq*lambda_lambda1;
    if sum(eig(sigaq)<0)>0
        alpha_lambda = alpha_lambda_old;
        muaq         = muaq_old;
        sigaq        = sigaq_old;
    end
    if sum(eig(siglq)<0)>0
       lambda_lambda = lambda_lambda_old;
       mulq          = mulq_old;
       siglq         = siglq_old;
    end
    
    as = sigma_lambda(1);
    bs = sigma_lambda(2);
    ag = gamma_lambda(1);
    bg = gamma_lambda(2);
    if as < 0
       as = as_old; 
    end
    if bs < 0
       bs = bs_old;
    end
    if ag < 0
       ag = ag_old; 
    end
    if bg < 0
       bg = bg_old; 
    end
    sigma_lambda = [as;bs];
    gamma_lambda = [ag;bg];
    lb = lb/S;
    lbc = [lbc,lb];
    disp(['lb: ',num2str(lb)]);
end



