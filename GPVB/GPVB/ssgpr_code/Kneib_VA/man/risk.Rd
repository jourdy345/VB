\name{risk}
\alias{risk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Risk for mean and quantile estimations
}
\description{
This functions can be used to evaluate the estimation of mean or quantile regression models. It or returns the MSE or the risk function calculated by the checkfunction.

}
\usage{
risk(y, pred, quant = FALSE, tau = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
  numeric vector of data  
%%     ~~Describe \code{y} here~~
}
  \item{pred}{
  numeric value or vector of the prediction for y
%%     ~~Describe \code{pred} here~~
}
  \item{quant}{
  logical. If TRUE the checkfunction is calculated (the default is FALSE)
%%     ~~Describe \code{quant} here~~
}
  \item{tau}{
  numerical. The quantile for which the risk is to be calculated. This value has to be between 0 and 1 (the default value is 0.5).
%%     ~~Describe \code{tau} here~~
}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
y <- rnorm(1000,3,5)
x <- mean(y)
## MSE of the mean
risk(y,x)

x <- quantile(y, .4)
###
risk(y,x,quant=TRUE, tau=.4)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ risk }
%\keyword{ mse }% __ONLY ONE__ keyword per line
