function [llh,mu_n,var_n] = bslogl_est(p,syn_n,S_i,S_yobs)
    %%% Requires n summary statistics in S_y
   
    mu_n = mean(S_i,1);
    
    
    var_n = zeros(size(S_i,2),size(S_i,2));
    for loop = 1:size(S_i,1)
        var_n = var_n + 1/(syn_n-1)* (S_i(loop,:) - mu_n)'*(S_i(loop,:) - mu_n);
    end

    
    psi_syn = 0;
    for loop = 1:p
        psi_syn = psi_syn + psi((syn_n - loop)/2); 
    end
   % llh = -p/2*log(2*pi) - 0.5*logdet(var_n) - 0.5*p*log((syn_n - 1)/2) + 0.5*psi_syn - 0.5*(syn_n-p-2)/(syn_n-1)*(S_yobs - mu_n)*inv(var_n)*(S_yobs - mu_n)' + 0.5*p/(syn_n-1) ;
    
     llh = -p/2*log(2*pi) - 0.5*logdet(var_n) - 0.5*p*log((syn_n - 1)/2) + 0.5*psi_syn   - 0.5*(syn_n-p-2)/(syn_n-1)*(S_yobs - mu_n)/(var_n)*(S_yobs - mu_n)' + 0.5*p/(syn_n);
end