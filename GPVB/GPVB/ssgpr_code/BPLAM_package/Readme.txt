We have three files include codes and a simulation example. In the code files, there are main codes used in the paper. In the simulation file, there is a simple simulation example.

############################################################################################

Code files:

1. VB_MC.R:       main function (BPLAM_VB) for the VB_MC algortihm
2. VB_Laplace.R:  main function (BPLAM_VB) for the VB_Laplace algortihm

Usage: 
1. BPLAM_VB(X,Y,l) 
2. BPLAM_VB(X,Y,l) 

Input (same for two functions): 

2. X is the covariate matrix, which is a n*p matrix with n the sample size and p the dimension of predictors.
3. Y is the response, which is a n*1 matrix.
4. l is the length of each subinterval, which can partition [0,1] in to k=length(seq(0, 1, l))-2 subintervals. We use 0.2, resulting in a total of K=6 bases

Output (same for two functions):

mu_u/xi_u: estimated mean and variance of \mu (scalars).

mu_a/xi_a: estimated mean and variance of \alpha (vectors of length p).

mu_b/xi_b: estimated mean and covariance matrix of \beta (matrices of size K*p / K*(K*p))

p_gamma.a: estimated probabilities \gamma^(alpha)_j=1, j=1...p (a vector of length p).

p_gamma.b: estimated probabilities \gamma^(beta)_j=1, j=1...p (a vector of length p).

A_sigma.sq/B_sigma.sq: estimated shape/scale parameters of \sigma^2_j, j=1...p (vectors of length p).

A_tau.sq/B_tau.sq: estimated shape/scale parameters of \tau^2_j, j=1...p (vectors of length p).

A_delta.sq/B_delta.sq: estimated shape/scale parameters of \delta^2 (scalars).


-----------------------------------------------------------------------------------------------------------

Simulation file: 

1. simulation.R

With the code in "simulation.R", data set described in the simulation example 1 can be generated and fitted (before using this code, please put R files "VB_MC.R" and "VB_Laplace.R" into your working directory and source either one accordingly). 
