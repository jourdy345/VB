load('/Users/daeyounglim/Desktop/Github/611_graphical/daeyoung/SSUR_EMVS/MATLAB/SSUREMVS/MATLAB_MCMC/Data/Simu_p6.mat') % Load data used in the paper; Other variables contained in 
% Simu_p6.mat are the results based on the following SSUR model analysis on
% the data.

[n,p,T] = size(F); 

n1 = 6;     % True regression coefficients
beta1=[1.3 0  -0.5  zeros(1,n1)];
beta2=[0.9 -.3  0.5  zeros(1,n1)];
beta3=[1  0.5  0.7 zeros(1,n1)];
beta4=[0.8 -0.6 zeros(1,n1)];
beta5=[1 0.7,zeros(1,n1)];
beta6=[1.1  0.6,zeros(1,n1)];
beta_true = [beta1'; beta2'; beta3'; beta4';beta5'; beta6'];
z1_true = find(abs(beta_true)>.001);  % z1_true: true subset of variables;


phi = .6;  % True error covariance matrix
Cov_true = 1/(1-phi^2).* toeplitz(phi.^([0:p-1]));   
Lambda=inv(Cov_true); adj_true = abs(Lambda)>.001; % adj_true: true adjacency matrix;


% initialize variational parameters
V = wishrnd(eye(p),4);
k = 4;
