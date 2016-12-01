% log posterior
lp = @(x)(x-log(1+exp(x)));
xp=-5:0.01:5;
plot(xp,lp(xp))

% define approximaton
m = 0
v = 2
lq = @(x)(-0.5*log(2*pi*v) - (0.5/v)*(x-m).^2);

% control variates stuff
dq1 = @(x)bsxfun(@minus,x,m);
dq2 = @(x)bsxfun(@minus,-0.5*x.^2,-0.5*(m^2+v));
dq = @(x)([dq1(x), dq2(x)]);
h1 = @(x)(((dq1(x)-mean(dq1(x)))'*bsxfun(@minus,dq(x),mean(dq(x))))/(length(x)-1) - [v, -m*v]);
h2 = @(x)(((dq2(x)-mean(dq2(x)))'*bsxfun(@minus,dq(x),mean(dq(x))))/(length(x)-1) - [-m*v, (m^2+0.5*v)*v]);
z = @(x)((x-m)/sqrt(v));
h1g = @(x)(-v*x + m*v);
h2g = @(x)(bsxfun(@minus,-m*v*[ones(length(x),1), -x]-(0.5*v^1.5)*[z(x), -x.*z(x)],[-m*v, (m^2+0.5*v)*v]));
grad_q_min_p = @(x)((m/v)-(1/v)*x - (1-1./(1+exp(-x))));
kl_grad = @(x)[v*grad_q_min_p(x), -m*v*grad_q_min_p(x)-(0.5*v^1.5)*grad_q_min_p(x).*z(x)];


% pre-allocate containers
nrit = 1e5;
gradest = zeros(2,11,nrit);

% fix random number seed
rng(12345)


% do repeats
for i=1:nrit
    
    % samples
    xi = m + sqrt(v)*randn(50,1);
    xi1 = xi(1:25);
    xi2 = xi(26:end);
    
    % simple approximation
    dqi = dq(xi);
    dpq = lq(xi)-lp(xi);
    gradest(:,1,i) = dqi'*dpq/50;
    
    % covariance approximation
    gradest(:,2,i) = bsxfun(@minus,dqi,mean(dqi))'*(dpq-mean(dpq))/49;
    
    
    % control variates, first method
    h1vec = zeros(5,2);
    h2vec = zeros(5,2);
    fvec = zeros(5,2);
    for j=1:5
        xij = xi1((1+(j-1)*5):(j*5));
        h1vec(j,:) = h1(xij);
        h2vec(j,:) = h2(xij);
        dpqij = lq(xij)-lp(xij);
        dqij = bsxfun(@minus,dq(xij),mean(dq(xij)));
        fvec(j,:) = (dpqij-mean(dpqij))'*dqij/4;
    end
    h1vec = bsxfun(@minus,h1vec,mean(h1vec));
    h2vec = bsxfun(@minus,h2vec,mean(h2vec));
    fvec = bsxfun(@minus,fvec,mean(fvec));
    alpha_1 = (h1vec'*h1vec)\(h1vec'*fvec(:,1));
    alpha_2 = (h2vec'*h2vec)\(h2vec'*fvec(:,2));
    cv1 = h1(xi2);
    cv2 = h2(xi2);
    dqi2 = dq(xi2);
    ddqi2 = bsxfun(@minus,dqi2,mean(dqi2));
    dpq2 = lq(xi2)-lp(xi2);
    ddpq2 = dpq2 - mean(dpq2);
    gradest(:,3,i) = ddqi2'*ddpq2/24 - [cv1*alpha_1; cv2*alpha_2];
    
    
    % control variate coefficients, 'stochastic linear regression'
    dqi1 = dq(xi1);
    ddqi1 = bsxfun(@minus,dqi1,mean(dqi1));
    dpq1 = lq(xi1)-lp(xi1);
    ddpq1 = dpq1 - mean(dpq1);
    alpha_lr = (ddqi1'*ddqi1)\(ddqi1'*ddpq1);
    gradest(:,4,i) = ddqi2'*ddpq2/24 - [cv1; cv2]*alpha_lr;
    
    
    % control variates, Ranganath et al.
    dq1i1 = dq1(xi1);
    ddq1i1 = dq1i1 - mean(dq1i1);
    dq2i1 = dq2(xi1);
    ddq2i1 = dq2i1 - mean(dq2i1);
    f1i1 = dq1i1.*dpq1;
    df1i1 = f1i1 - mean(f1i1);
    f2i1 = dq2i1.*dpq1;
    df2i1 = f2i1 - mean(f2i1);
    a1 = ddq1i1'*df1i1/(ddq1i1'*ddq1i1);
    a2 = ddq2i1'*df2i1/(ddq2i1'*ddq2i1);
    dq1i2 = dq1(xi2);
    dq2i2 = dq2(xi2);
    f1i2 = dq1i2.*dpq2;
    f2i2 = dq2i2.*dpq2;
    gradest(:,5,i) = [mean(f1i2-a1*dq1i2); mean(f2i2-a2*dq2i2)];
    
    
    % 2nd order taylor of log p around m, Paisley et al.
    lp_at_m = lp(m);
    prob_at_m = 1./(1+exp(-m));
    taylor_prec = prob_at_m*(1-prob_at_m);
    taylor_meanprec = m*taylor_prec + (1-prob_at_m);
    taylor_approx = @(x)(lp_at_m + taylor_meanprec*(x-m) - 0.5*taylor_prec*(x-m).^2);
    taylor_cv = @(x)bsxfun(@times,dq(x),taylor_approx(x));
    t_cv_1 = taylor_cv(xi1);
    t_cv_1 = bsxfun(@minus,t_cv_1,mean(t_cv_1));
    a = sum(sum(t_cv_1.*bsxfun(@times,dq(xi1),lp(xi1))))/sum(sum(t_cv_1.^2));
    paisley_f_2 = lp(xi2)-a*taylor_approx(xi2);
    gradest(:,6,i) = [v, -m*v; -m*v, (m^2+0.5*v)*v]*[(m/v)-a*(taylor_meanprec+taylor_prec*m); (1/v)-a*taylor_prec] ...
        - [mean(dq1i2.*paisley_f_2); mean(dq2i2.*paisley_f_2)];
    
    
    % gradient, Kingma & Welling
    gradest(:,7,i) = mean(kl_grad(xi))';
    
    
    % biased regression estimate 1
    ddqi = bsxfun(@minus,dqi,mean(dqi));
    alpha_lr = (ddqi'*ddqi)\(ddqi'*dpq);
    gradest(:,8,i) = [v, -m*v; -m*v, (m^2+0.5*v)*v]*alpha_lr;
    
    
    % biased regression estimate 2
    gt = [1, -mean(xi)];
    covmat = [v*gt; ...
        -m*v*gt-(0.5*v^1.5)*mean([z(xi) -xi.*z(xi)])];
    alpha_lr = covmat\gradest(:,7,i);
    gradest(:,9,i) = [v, -m*v; -m*v, (m^2+0.5*v)*v]*alpha_lr;
    
    
    % biased regression estimate 3
    prob_at_x = 1./(1+exp(-xi));
    grad_at_x = 1-prob_at_x;
    hess_at_x = -prob_at_x.*grad_at_x;
    alpha_lr = [(m/v); 1/v] - [mean(grad_at_x)-mean(hess_at_x)*mean(xi); -mean(hess_at_x)];
    gradest(:,10,i) = [v, -m*v; -m*v, (m^2+0.5*v)*v]*alpha_lr;
    
    
    
    % control variates, first method + gradient
    h1vec = h1g(xi1);
    h2vec = h2g(xi1);
    fvec = kl_grad(xi1);
    h1vec = bsxfun(@minus,h1vec,mean(h1vec));
    h2vec = bsxfun(@minus,h2vec,mean(h2vec));
    fvec = bsxfun(@minus,fvec,mean(fvec));
    alpha_1 = (h1vec'*h1vec)\(h1vec'*fvec(:,1));
    alpha_2 = (h2vec'*h2vec)\(h2vec'*fvec(:,2));
    cv1 = h1g(xi2);
    cv2 = h2g(xi2);
    gradest(:,11,i) = mean(kl_grad(xi2) - [cv1*alpha_1, cv2*alpha_2])';
    
end


% get mean gradients
meangrads = mean(gradest,3)
vargrads = sum(sum(bsxfun(@minus,gradest,meangrads).^2,3),1)/(nrit-1)

% biases
sq_bias_my_grad1 = sum((meangrads(:,8)-meangrads(:,3)).^2)
sq_bias_my_grad2 = sum((meangrads(:,9)-meangrads(:,3)).^2)
sq_bias_my_grad3 = sum((meangrads(:,10)-meangrads(:,3)).^2)