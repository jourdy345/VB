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
invsig_q = inv(sig_q);
T_mat = chol(invsig_q,'lower');
T_mat_transpose = T_mat';
T_mat_prime = T_mat;
for k = 1:(n+3)
   T_mat_prime(k,k) = log(T_mat(k,k));
end
niter = 10000;
count = 0;
lbc = [];
mu_q_storage = zeros(n+3,niter);

% start inference
for i = 1:niter
   count = count+1;
   rho = 1/(5+count);
   T_mat_old = T_mat;
   T_mat_transpose_old = T_mat_transpose;
   T_mat_prime_old = T_mat_prime;
   mu_q_old = mu_q;
   
   s = randn(n+3,1);
   solveTs = T_mat_transpose\s;
   theta = mu_q+solveTs;
   b = theta(1:n);
   alpha = theta(n+1);
   lam = theta(n+2);
   psi = theta(n+3);
   phi = exp(psi)/(1+exp(psi));
   sig = exp(alpha);
   Ts = T_mat*s;
   g_b1 = -(1-phi*phi)*b(1)+phi*(b(2)-phi*b(1))-sig/2+sig/2*y(1)*y(1)*exp(-lam-sig*b(1));
   g_b = phi*(b(3:n)-phi*b(2:(n-1)))-(b(2:(n-1))-phi*b(1:(n-2)))-sig/2+sig/2*y(2:(n-1)).*y(2:(n-1)).*exp(-lam-sig*b(2:(n-1)));
   g_bn = -(b(n)-phi*b(n-1))-sig/2+sig/2*y(n)*y(n)*exp(-lam-sig*b(n));
   g_a = 0.5*sum(y.*y.*b.*exp(alpha-lam-sig*b))-sig/2*sum(b)-alpha/sig_a2;
   g_l = -n/2+0.5*sum(y.*y.*exp(-lam-sig*b))-lam/sig_l2;
   g_p = exp(psi)/((1+exp(psi))*(1+exp(psi)))*(phi*b(1)*b(1)-phi/(1-phi*phi)+sum((b(2:n)-phi*b(1:(n-1))).*b(1:(n-1))))-psi/sig_p2;
   g_mu = [g_b1;g_b;g_bn;g_a;g_l;g_p]+Ts;
   
   % update mu
   mu_q = mu_q+rho*g_mu;
   g_T = -solveTs*(T_mat\g_mu)';
   for k = 1:(n+3)
      g_T(k,k) = g_T(k,k)*T_mat(k,k); 
   end
   T_mat_prime = T_mat_prime+rho*g_T;
   T_mat = T_mat_prime;
   for k = 1:(n+3)
      T_mat(k,k) = exp(T_mat_prime(k,k));
   end
   T_mat_transpose = T_mat';
   T_mat = tril(T_mat);
   sig_q = T_mat*T_mat';
   if sum(sum(~isfinite(sig_q)))>0
       T_mat = T_mat_old;
       T_mat_transpose = T_mat_transpose_old;
       T_mat_prime = T_mat_prime_old;
       mu_q = mu_q_old;
   elseif any(eig(sig_q)<0)
       T_mat = T_mat_old;
       T_mat_transpose = T_mat_transpose_old;
       T_mat_prime = T_mat_prime_old;
       mu_q = mu_q_old;
%    end
%        T_mat = T_mat_old;
%        T_mat_transpose = T_mat_transpose_old;
%        T_mat_prime = T_mat_prime_old;
%        mu_q = mu_q_old;
   end
   lb = logh(n,alpha,b,y,lam,phi,psi,sig_a2,sig_l2,sig_p2)+(n+3)/2*log(2*pi)-logdet(T_mat)+0.5*sum(s.*s);
   lbc = [lbc,lb];
   mu_q_storage(:,i) = mu_q;
   disp(['count: ',num2str(count),', lb: ',num2str(lb)]);
end


