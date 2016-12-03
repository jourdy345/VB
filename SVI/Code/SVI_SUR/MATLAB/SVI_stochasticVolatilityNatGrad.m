n = 200;
phi_true = 0.3;
lam_true = 1.24;
sig_true = -1.73;
b1_true = randn(1,1)*sqrt(1/(1-phi_true*phi_true));
b_true = zeros(n,1);
b_true(1) = b1_true;
y = zeros(n,1);
y(1) = randn(1,1)*exp(0.5*(lam_true+sig_true*b_true(1)));
for t = 2:n
   b_true(t) = randn(1,1)+phi_true*b_true(t-1);
   y(t) = randn(1,1)*exp(0.5*(lam_true+sig_true*b_true(t)));
end

sig_a2 = 1;
sig_l2 = 1;
sig_p2 = 1;

mu_q  = zeros(n+3,1);
sig_q = eye(n+3);
D_mat = duplication_matrix_full(n+3);
D_plus = duplication_matrix_MP_inv(n+3);
invsig_q = invpd(sig_q);
gauss_lambda1 = sig_q\mu_q;
gauss_lambda2 = -0.5*D_mat'*invsig_q(:);
gauss_lambda = [gauss_lambda1;gauss_lambda2];
% T_mat = chol(invsig_q,'lower');
% T_mat_transpose = T_mat';
% T_mat_prime = T_mat;
% for k = 1:(n+3)
%    T_mat_prime(k,k) = log(T_mat(k,k));
% end
niter = 10000;
S = 200; % number of simulations
count = 0;
lbc = [];
mu_q_storage = zeros(n+3,niter);

% start inference
for i = 1:niter
    gauss_lambda_old = gauss_lambda;
    mu_q_old = mu_q;
    sig_q_old = sig_q;
    rho = 1/(5+i);
    M_mat = 2*D_plus*kron(mu_q,eye(n+3));
    S_mat = 2*D_plus*kron(sig_q,sig_q)*D_plus';
    aux = S_mat\M_mat;
    gauss_invfisherinformation = [invsig_q+M_mat'*aux -aux';-aux eye(size(S_mat))/S_mat];
    theta = mvnrnd(mu_q,sig_q,S);
   
end