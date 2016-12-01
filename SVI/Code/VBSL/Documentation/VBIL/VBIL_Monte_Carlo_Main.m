%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Normal model example 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This code uses ABC likelihood
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generating the data
clear all
rng(100);
theta_true = normrnd(0,1,1,1);
n = 4;
sum_length = n;


y_obs = zeros(n,1);
S_obs = y_obs;

eps = 0.1*(n+1)/(n-0.1); 
Sign = eps*eye(sum_length);
Siginv_ABC = eye(sum_length)/Sign; %%% this is for the ABC likelihood, not the inverse of the Sig used
N_start = 50;
sig2_opt = 0.1;

Monte = 1;  %%% this is just an indicator when saving the file

%% Initialization

d = 1;
D_d = duplication_matrix_full(d);
D_plus = duplication_matrix_MP_inv(d);

tic
S = 100;
theta = mean(y_obs);
mu = theta; 
Sig = eye(d);
aux_matrix = eye(d)/Sig;
lambda1 = Sig\mu;
lambda2 = -1/2*D_d'*aux_matrix(:);
lambda = [lambda1;lambda2];

%% Step 1 of Algorithm 1 to estimate c

Sig = -1/2*eye(d)/reshape(D_plus'*lambda2,d,d);
mu = Sig*lambda1;
gra_log_q_lambda = zeros(S,d+d*(d+1)/2); 
g_lambda = zeros(S,d+d*(d+1)/2); 
LB = 0;
log_llh = zeros(S,1);
scalar = zeros(S,1);
var_z = zeros(S,1);
N = zeros(S,1);


rqmc = gen_Sobol(ceil(log2(S)),d)'; % generate randomised QMC numbers
parfor s = 1:S  
    
    U_normal = norminv(rqmc(s,:))'; 
    theta = mu+chol(Sig,'lower')*U_normal;
 
    
    [llh var_z(s) N(s)] = log_ABC_est(theta,Siginv_ABC,n,S_obs,sig2_opt,N_start);
    log_llh(s) = llh;
    
    prior = log(mvnpdf(theta,0,1*eye(d)));
    logqtheta = log(mvnpdf(theta,mu,Sig));

    scalar(s) = llh+prior - logqtheta;
    
    gra_log_q_lambda(s,:) = [theta-mu;vech(theta*theta'-Sig-mu*mu')]';
    g_lambda(s,:) = gra_log_q_lambda(s,:)*scalar(s);    
    

end
c12 = zeros(1,d+d*(d+1)/2); 
parfor i = 1:d+d*(d+1)/2
    aa = cov(g_lambda(:,i),gra_log_q_lambda(:,i));
    c12(i) = aa(1,2)/aa(2,2);
end

%% Data storage

niter = 100;iter = 1;
countdown = 1;
stop = false;
stop_converge = false;



storemean = zeros(niter+countdown,4);
storevar = zeros(niter+countdown,4);

storemean(iter,:) = mu;
storevar(iter,:) = diag(Sig);

Nstore = zeros(niter+countdown,1);

S_check1 = 20;
%S_check2 = 20;
grad_sd = zeros(S_check1,d+d*(d+1)/2,niter+countdown);



%% Iterate main algorithm until convergence
while ~stop   
    
    lambda_old = lambda;
    g_lambda_for_c = zeros(S,d+d*(d+1)/2);
    log_llh = zeros(S,1);
    
    
    
    %% Checking standard deviaton of gradient estimate %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
% 
%     gra_log_q_lambda_check = zeros(S,d+d*(d+1)/2); 
%     g_lambda_check = zeros(S,d+d*(d+1)/2); 
%  
%     for ss = 1:S_check1   
%     rqmc = gen_Sobol(ceil(log2(S)),d)';
%     parfor s = 1:S
%         
%         U_normal = norminv(rqmc(s,:))'; 
%         theta = mu+chol(Sig,'lower')*U_normal;
%         
%         prior = log(mvnpdf(theta,0,1*eye(d)));
% 
%         [llh, var_z_check, N_check] = log_ABC_est(theta,Siginv_ABC,n,S_obs,sig2_opt,N_start);
%                
%         %logqtheta = 55; 
%         logqtheta = log(mvnpdf(theta,mu,Sig));
% 
%          
%         scalar_check = llh+prior - logqtheta;
%         gra_log_q_lambda_check = [theta-mu;vech(theta*theta'-Sig-mu*mu')]';
%         g_lambda_check(s,:) = gra_log_q_lambda_check.*(scalar_check-c12);    
%       
%     end
%         Y_check = mean(g_lambda_check)';
%         
%         mat_M = 2*D_plus*kron(mu,eye(d));
%         mat_S = 2*D_plus*kron(Sig,Sig)*D_plus';
%         aux = mat_S\mat_M;
%         I_Fisher_lambda_inv = [eye(d)/Sig+mat_M'*aux -aux';-aux eye(size(mat_S))/mat_S];
%         grad_sd(ss,:,iter) = I_Fisher_lambda_inv*Y_check;
%     end    
    
    %% Main Code - See Step 2 Algorithm 1 %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    rqmc = gen_Sobol(ceil(log2(S)),d)';
    parfor s = 1:S    
        
        U_normal = norminv(rqmc(s,:))'; 
        theta = mu+chol(Sig,'lower')*U_normal;


        %% Calculation of VBIL 
        
        [llh var_z(s) N(s)] = log_ABC_est(theta,Siginv_ABC,n,S_obs,sig2_opt,N_start);
        log_llh(s) = llh;

        %% Calculation of gradient estimate
        
        prior = log(mvnpdf(theta,0,1*eye(d)));
        logqtheta = log(mvnpdf(theta,mu,Sig));
        
        %%% h(theta,z)
        scalar(s) = llh+prior - logqtheta;
        gra_log_q_lambda(s,:) = [theta-mu;vech(theta*theta'-Sig-mu*mu')]';
        g_lambda(s,:) = gra_log_q_lambda(s,:).*(scalar(s)-c12);    
        g_lambda_for_c(s,:) = gra_log_q_lambda(s,:)*scalar(s);  
    end
    c12 = zeros(1,d+d*(d+1)/2); 
    parfor i = 1:d+d*(d+1)/2
        aa = cov(g_lambda_for_c(:,i),gra_log_q_lambda(:,i));
        c12(i) = aa(1,2)/aa(2,2);
    end
    Y = mean(g_lambda)';
    
    mat_M = 2*D_plus*kron(mu,eye(d));
    mat_S = 2*D_plus*kron(Sig,Sig)*D_plus';
    aux = mat_S\mat_M;
    I_Fisher_lambda_inv = [eye(d)/Sig+mat_M'*aux -aux';-aux eye(size(mat_S))/mat_S];
    
    stepsize = 1/(5+iter/5);
    %lambda = (1-stepsize)*lambda+stepsize*(I_Fisher_lambda_inv*Y);
    lambda = lambda+stepsize*(I_Fisher_lambda_inv*Y);
    
    %%% Checking if Sig is positive semi definite
    
    lambda2 = lambda(d+1:d+d*(d+1)/2);           
    Sig = -1/2*eye(d)/reshape(D_plus'*lambda2,d,d);
    if sum(eig(Sig)<0)>0 
        lambda = lambda_old;  
    end      
    
    lambda1 = lambda(1:d);
    lambda2 = lambda(d+1:d+d*(d+1)/2);                  
    Sig = -1/2*eye(d)/reshape(D_plus'*lambda2,d,d);
    mu = Sig*lambda1;  %%% Is mu supposed to be updated?
    
    

    
    %LB(iter) = log(det(Sig))/2+mean(log_llh);
  %  LB(iter) = log(det(Sig))/2 + mean(scalar) + d/2*(1+ log(2*pi));
   LB(iter) = mean(scalar) ;
    LB(iter) = LB(iter)/n;
    lower_bound = LB'
    
      if iter>5
     %  if (abs(mean(LB(iter-3:iter))-mean(LB(iter-4:iter-1)))<0.0001||iter>niter) 
       if (iter>niter) 
           stop_converge = true; 
       end
  %     if (iter>niter) stop = true; end
    end
    
    if stop_converge == 1
       countdown = countdown - 1;
       if countdown == 0
          stop = true; 
       end
    end
    
    Nstore(iter) = mean(N);
    iter = iter+1;
    storemean(iter,:) = mu;
    storevar(iter,:) = diag(Sig);
    iter
    mu
    Sig
    mean(N)
    
end
CPU_VBIL = toc/60


color = 'r'

i = 1;
x = -2:.001:2;
yy_VBIL = normpdf(x,mu(i),sqrt(Sig(i,i)));
subplot(1,2,1)
plot(x,yy_VBIL,color,'LineWidth',2);
hold on;

subplot(1,2,2)
plot(LB,color,'LineWidth',2);
hold on;

subplot(1,2,1)
x = [-2:.01:2];
norm = normpdf(x,(n/(n+1))*mean(y_obs),sqrt(1/(1+n)));
plot(x,norm,'black','LineWidth',2)



