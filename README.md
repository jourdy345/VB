# Variational Bayes

You have visited my repository for variational inference! I am currently working on stochastic variational inference(SVI). The SVI codes are toy examples of Bayesian linear regression, one of which is stochastic gradient with no variance reduction and the other with variance reduction via control variate. There is a vast literature on variational inference out there so grab one paper and try to read it through. You can use my codes as examples of those models. Feel free to ask me if you have any question. I will try my best to reply.

<!-- Welcome to my repository for variational inference! The repository you've just stumbled upon is not for public use/exhibit; it has been created entirely for personal purpose. I don't mind people using my code or reading my summary notes but please not that there is no manual. This README file has also been written to inform MYSELF and not to forget. So you're seeing my notebook basically.

### GP mixed C++ code compiling & linking
```bash
cd GP_mixed/Code
g++ -o GP main.cpp gaussianCDF.cpp blockDiag.cpp lowerBound.cpp eOfZ.cpp eOfZTZ.cpp logH.cpp sparseGPVBProbit.cpp -llapack -lblas -larmadillolarmadillo
```
Don't forget to link the added libraries.

## Simulation Code
There are three simulation codes:  

- `sim_GP.R`: This is the simulation code for the original code given by the author of __Variational inference for sparse spectrum Gaussian process regrsion, Linda S. L. Tan and David J. Nott__.
- `GPnormal_simulation.R`: Original model + random effects.
- `GPprobit_simulation.R`: Extension to probit model with random effects.

1. `sim_GP.R`:  Function `sim_GP` takes in two parameters: `FUN` and `m`. `FUN` is a function the researcher defined and wants to test. `m` is the number of bases (sine/cosine).
2.  `GPnormal_simulation.R`: `sim_GPnormal` takes at most 4 parameters: `FUN`, `m`, `intercept`, and `draw`. `intercept` and `draw` are Boolean variables. Their default values are `TRUE`. `intercept` can be set to `FALSE` if the researcher does not want to include the intercept. `draw` is set to `FALSE` if the researcher does not want the function to generate the resulting plot at the end.
3. `GPprobit_simulation.R`: `sim_GPprobit` takes at most three parameters: `FUN`, `m`, and `intercept`. Same as above.

####Example

```R
	my_fun = function(x) 2 * x - 1
	fit = sim_GP(my_fun, 10)
```
 
```R
	my_fun = function(x) tanh(4*x-2)
	fit = sim_GPnormal(my_fun, 10)
```

```R
	my_fun = function(x) 0.5*exp(x*2)-2
	fit = sim_GPprobit(my_fun, 10)
```

### Note
----
If a user-defined function is provided, `sim_GPnormal` and `sim_GPprobit` both generate data according to the model `f(x) + Zu + ε` where `ε ~ N(0,2)` and `u ~ N(0,3)`. `f(x)` is approximated by `Σ( a*sin(2πsx) + b*cos(2πsx) )`. -->