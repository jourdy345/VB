\name{lin}
\alias{lin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Linear effects for Bayesian inference
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Simplifies the access to prior parameters given in a formula type object.
}
\usage{
lin(x, mu0=0, sig0=10000)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
 \item{x}{
 data
 }
 \item{mu0}{
 mean for the gaussian prior for the linear effect
 }
 \item{sig0}{
 variance for the gaussian prior for the linear effect
 }
 }


\value{
The vector returned only comprises the data. The attributes are the class of the effect("linear"), and the prior parameters.  
}

\author{
Elisabeth Waldmann
%%  ~~who you are~~
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{bsplines}, \link{spatial}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
x <- runif(100,-3,3)
Z <- lin(x)
attr(Z, "class")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
