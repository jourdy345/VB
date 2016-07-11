#########################################################################
### This is a sample code for demonstration that works on my computer ###
###         Change the path variable and file name variables          ###
########################################################################3
path           <- '~/Desktop/VB/GPVB/SSGP/Miguel_datasets/' # make sure it ends with a forward slash
fileName_X_tr  <- 'X_tr.txt' # don't leave out the file extension
fileName_T_tr  <- 'T_tr.txt'
fileName_X_tst <- 'X_tst.txt'
fileName_T_tst <- 'T_tst.txt'
fit            <- compareSSGPvsBSAR(data = 'pendulum', fit = 'test', fileName_X_tr = fileName_X_tr, fileName_T_tr = fileName_T_tr, fileName_X_tst = fileName_X_tst, fileName_T_tst = fileName_T_tst)

## Or choose to simulate a function
fit_simulate   <- compareSSGPvsBSAR(data = 'simul') # R will ask for a function definition. Use proper R syntax.
