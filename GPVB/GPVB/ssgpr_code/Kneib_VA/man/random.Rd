 \name{random}
\alias{random}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Random effects
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Generates a basis matrix for a random effect facilitating the access to parameters in a formula type object
for Bayesian inference (such as the prior parameters)
}

\usage{
random(ind, a_1=1, b_1=1e-5)
}
%- maybe also 'usage' for other objects documented here.

  \arguments{
 \item{ind}{
 information about the individuals or group}
 \item{a_1}{
 shape parameter for the gamma-prior for the smoothing parameter}
 \item{b_1}{
 scale parameter for the gamma-prior for the smoothing parameter}
 }


\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
The matrix returned is the basis matrix. The attributes are
  K (the neighbouring matrix) and the class of the effect ("random") as well as the prior parameters a_1 and b_1.
  %% ...
}
\references{
%% ~put references to the literature/web site here ~
Fahrmeir L., Kneib T. and Lang S. (2009) \emph{Regression} Springer, New York
}
\author{
Elisabeth Waldmann

\email{ewaldma@uni-goettingen.de}

}

\seealso{
\link{lin}, \link{bsplines}, \link{spatial}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
 