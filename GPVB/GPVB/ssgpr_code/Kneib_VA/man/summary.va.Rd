 \name{summary.va}
\alias{summary.va}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary for VA-Regression
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Summarizes the results of the VA-Regerssion
}

\usage{
\method{summary}{va}(object, ...)
}
%- maybe also 'usage' for other objects documented here.

  \arguments{
 \item{object}{
 object of va-type}
 \item{...}{Additional arguments}
 }




\details{
This function returns the formula-object, the quantile the regression was done for (if \code{quant=TRUE}),
a table of the linear coefficients, containing the standard deviation as well as the 
0.025 and the 0.975-quantile of the corresponding Gaussian distribution. 
The nonlinear, spatial and random effects are mentioned with their meaning.
}


\author{
Elisabeth Waldmann \cr
\email{ewaldma@uni-goettingen.de}

}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
 