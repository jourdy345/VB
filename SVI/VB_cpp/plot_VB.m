x = csvread('x_f1');
y = csvread('y_f1');
fx = csvread('fx_f1');
L_post = csvread('L_post_f1');
thetac = csvread('thetaContainer_f1');


vphi = sqrt(2)*cos(x*(pi*(1:35)));

Sig = L_post * L_post';
Sig = Sig(1:35,1:35);
theta_est = thetac(1:35,10000);
fx_est = vphi*theta_est;

n = 200;
CItheta = mvnrnd(theta_est',Sig,5000);
Fullfx = zeros(n,5000);
for i = 1:5000
   fx_ci = vphi(:,1:length(CItheta(i,:)))*CItheta(i,:)';
   Fullfx(:,i) = fx_ci;
end
CIfx = quantile(Fullfx,[0.025,0.975],2);

lower = CIfx(:,1);
upper = CIfx(:,2);

[~,o] = sort(x);

h = ciplot(lower(o),upper(o),x(o),'k');
hold on
set(h,'facealpha',.3);
plot(x(o),fx(o),x(o),fx_est(o));
legend('','true','estimate');
hold off


for i = 1:39
   subplot(6,7,i)
   plot(thetac(i,:));
end