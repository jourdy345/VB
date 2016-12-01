function [llh,var_llh,N] = log_ABC_est(theta,Siginv_ABC,n,S_obs,sig2_opt,N_start)
%N_start = 50;
%N_inc = 50;
N_inc = 50;
N = N_start;
lw = zeros(N,1);
for i = 1:N
    
    y = normrnd(theta,1,n,1);
      
    S = y;
    %S = mean(y);
    
    lw(i) = -1/2*(S-S_obs)'*Siginv_ABC*(S-S_obs) -length(S)/2*log(2*pi) - 0.5*log(det(inv(Siginv_ABC)));
end
max_lw = max(lw);    
var_llh = sum(exp(2*(lw-max_lw)))/(sum(exp(lw-max_lw)))^2-1/N;    
 while (var_llh>sig2_opt)&&(N<30000)
     lw_inc = zeros(N_inc,1);
     for i = 1:N_inc
         y = normrnd(theta,1,n,1);
         S = y;
         lw_inc(i) =  -1/2*(S-S_obs)'*Siginv_ABC*(S-S_obs) -length(S)/2*log(2*pi) - 0.5*log(det(inv(Siginv_ABC)));
     end
     lw = [lw;lw_inc];
     N = N+N_inc;
     max_lw = max(lw);    
     var_llh = sum(exp(2*(lw-max_lw)))/(sum(exp(lw-max_lw)))^2-1/N;    
 end

llh = max_lw+log(sum(exp(lw-max_lw)))-log(N);          %%% More stable than just log(sum(exp(lw)))-log(N)

end
