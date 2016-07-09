 \name{spatial}
\alias{spatial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gaussian random fields
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Generates a basis matrix for a Gaussian random field facilitating the access to parameters in a formula type object
for Bayesian inference (such as the prior parameters)
}

\usage{
spatial(geo, map,a_1=1, b_1=.0001)
}
%- maybe also 'usage' for other objects documented here.

  \arguments{
 \item{geo}{
 spatial information}
 \item{map}{
 the corresponding map (given as .gra or .bnd) }
 \item{a_1}{
 shape parameter for the gamma-prior for the smoothing parameter}
 \item{b_1}{
 scale parameter for the gamma-prior for the smoothing parameter}
 }



\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
The matrix returned is the basis matrix. The attributes are
  K (the neighbouring matrix) and the class of the effect ("spatial") as well as the prior parameters a_1 and b_1.
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{lin}, \link{bsplines}
}
\examples{
germany <- read.bnd(system.file("examples/germany.bnd", package="BayesX"))
gramap<-bnd2gra(germany)
centroids <- get.centroids(germany)
geo_ind <- c(sample(1:310,size=1000,replace=TRUE))
geo <- row.names(centroids)[geo_ind]
Z <- spatial(geo, gramap)
K <- attr(Z, "K")
attr(Z, "class")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
