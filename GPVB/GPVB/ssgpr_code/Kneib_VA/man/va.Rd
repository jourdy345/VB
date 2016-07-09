\name{va}
\alias{va}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of Parameters of Bayesian Additive (Quantile) Regression
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
function to estimate parameters of bayesian additive (quantile) regression
}
\usage{
va(formula, data = NULL, intercept = TRUE, quant = FALSE, ta = 0.5,
k = 100, sto = 1e-04,
a_0 = 1, b_0 = 1e-05, 
plotmse = FALSE, norm = TRUE, pred = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
  a formula object, with the response on the left of a ~ operator, and the terms, separated by + operators, on the right, by default calculated as linear effects. The effects can be changed by using different functions, \link{bsplines}, \link{spatial}, \link{lin} or \link{random}.
}


  \item{data}{
  a matrix containing the data
%%     ~~Describe \code{data} here~~
}
  \item{intercept}{
  logical value if an intercept should be included (the default value is TRUE)
%%     ~~Describe \code{intercept} here~~
}

  \item{quant}{
  logical value indicating if quantile regression is to be calculated (the default value is FALSE)
%%     ~~Describe \code{quant} here~~
}
  \item{ta}{
    numerical. The quantile to be estimated. This value has to be between 0 and 1 (the default value is 0.5).
%%     ~~Describe \code{ta} here~~
}
  \item{k}{
  maximal number of iterations used
%%     ~~Describe \code{k} here~~
}
  \item{sto}{
  stoping criterion, using the relative change of the risk to stop the program before reaching k iterations
%%     ~~Describe \code{sto} here~~
}
  \item{a_0}{
  shape-parameter for gamma-prior for modelprecision parameter
%%     ~~Describe \code{a_0} here~~
}
  \item{b_0}{
  rate-parameter for gamma-prior for modelprecision parameter
%%     ~~Describe \code{b_0} here~~
}
  \item{plotmse}{
  logical value, indicating if the risk should be plotted in each iteration (default value is FALSE). Note that the plotted risk refers to the data used for estimation (i.e. defaultwise the normalised data) and therefore might differ highly from the risk calculated afterwards.
}
  \item{norm}{
  logical value, indicating if the data are to be normalized before calculation (default value is TRUE)
}
  \item{pred}{
  matrix of data used for prediction. If there are no data for prediction, the prediction will be done on the data which were used for estimation
%%     ~~Describe \code{pred} here~~
}
  \item{\dots}{
Extra parameters needed
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use  
  \item{beta }{the estimated linear effects}
  \item{sig_bet}{the covariance matrix of the linear effects}
  \item{Xlinear}{a matrix containing all the linear information}
  \item{Z}{the designmatrix of the nonlinear effects}
  \item{gamma}{list of vectors of ALL estimated nonlinear effects (i.e. bsplines, spatial and random effects)}
  \item{sig_gam}{the covariance matrix of all nonlinear effects}
  \item{Xnonlinear}{a data matrix containing all the nonlinear information}
  \item{tau_inv}{a vector of the inverse of the smoothing parameters}
  \item{delta}{the model precision}
  \item{i}{the number of iterations until the stopping criterion or the maximal number of iterations was reached (note that there is a warning if k is reached before the stopping criterion is fulfilled)}
    \item{call}{formula given to the function}
  \item{prediction}{a vector with the predicted values for either the data the estimation was based on (if \code{pred=NULL}) or for the data given in pred}
  \item{tau}{the quantile the estimation was done for (of \code{quant = TRUE})}
  \item{residuals}{a vector of the difference between y and the prediction}
  \item{replace}{a vector which stores the numbers of the iterations in which the parameters in the distribution of the modelprecision were replaced in order to make the numerical integration possible. See section note for further information.}
  
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{

Ormerod, J.T. and Wand, M.P.
\emph{Explaining Variational Approximations}.
The American Statistician, (2010), 64, 140-153. 

Wand, M. P.; Ormerod, J. T.; Padoan, S. A.; and Fruhwirth, R., \emph{Variational Bayes for Elaborate Distributions},
\url{http://ro.uow.edu.au/cssmwp/56} (working paper).


}
\author{
Elisabeth Waldmann \cr
University of Goettingen \cr
\email{ewaldma@uni-goettingen.de}

}
\note{
If you estimate the coeffecients for a quantile regression, the approximated precision of the model is calculated via numerical integration using a subroutine in C. If values in this subroutine get too big, the step is repeated with the compensatory values for the parameters of the precision. This only happens for very sparse regions in the datasets (i.e. usually in the margins) and even in these cases the integral will normally be calculated for the real values before convergence. If the real values were not used throughout the whole iterating process, there will be a warning. If the results of your estimations are somehow akward even though there was no warning, you should check \code{va$replace}, which returns a vector, telling you in which iteration the parameters were replaced. If you set \code{plotmse=TRUE}, the points in which the parameters are replaced are marked in red. If they are replaced in an iteration close to the last iteration, you might want to think about the sparsity of the data in the quantile you tried to estimate (estimating the .01-quantile of a dataset containing 50 observations might not be very sensible anyway). If you want to use the estimations anyway, you should know that even if the point estimators (such as the beta-vector)  still seem trustworthy, the covariances get to big (you might consider taking the last estimation BEFORE the replacement to be the closest to the "truth").
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Graphical presentation can be done by using \link{plot.va}.
}
\examples{
#data(cars)
#va(dist~speed, data=cars)

## Example for using linear and nonlinear effects
set.seed(123)
x <- runif(100,-3,3)
x2 <- runif(100,-3,3)
x3 <- runif(100,-3,3)
x4 <- runif(100,0,3)
y <- 3 + 3*sin(x)+3*x2+4*x3+log(x4)+rnorm(100,0,1)

L<-va(y~bsplines(x)+x2 +x3 +bsplines(x4), quant=TRUE)
summary(L)
plot(L)

## Example using different quantiles in a heteroscedastic scenario

set.seed(123)
x <- runif(500,0,3)
y <- 3 + 1.5*x+rnorm(500,0,.3*x)

L1<-va(y~x, quant=TRUE, ta=.1)
L3<-va(y~x, quant=TRUE, ta=.3)
L5<-va(y~x, quant=TRUE, ta=.5)
L7<-va(y~x, quant=TRUE, ta=.7)
L9<-va(y~x, quant=TRUE, ta=.9)
plot(x,y, col="lightgrey")
lines(sort(x), L1$pred[order(x)])
lines(sort(x), L3$pred[order(x)])
lines(sort(x), L5$pred[order(x)])
lines(sort(x), L7$pred[order(x)])
lines(sort(x), L9$pred[order(x)])

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
