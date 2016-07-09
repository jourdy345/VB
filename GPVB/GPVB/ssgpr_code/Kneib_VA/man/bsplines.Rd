\name{bsplines}
\alias{bsplines}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
B-Splines Basis for Bayesian Inference
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Generates a b-spline basis matrix with attributes for Bayesian inference (such as the penalisation matrix and prior parameters). This function is based on the function \link{spline.des} from the package \pkg{splines}.
}
\usage{
bsplines(x, M = 20, degree = 3, order = 2, a_1=1, b_1=.00001, knots=NULL, mx=NULL)
}

  \arguments{
 \item{x}{data}
 \item{M}{number of knots}
 \item{degree}{the degree of splines}
 \item{order}{the order of the random walk prior}
 \item{a_1}{shape parameter for the gamma-prior for the smoothing parameter}
 \item{b_1}{scale parameter for the gamma-prior for the smoothing parameter}
 \item{knots}{knots at which the spline is to be evaluated}
 \item{mx}{value the spline is centered around, default is mean(x)}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
The matrix returned is the basis matrix. The attributes are
  K (the penalisation matrix with differences of order "order") and the class of the effect ("splines") as well as the prior parameters a_1 and b_1.
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
\link{lin}, \link{spatial}, \link{random}
}
\examples{
x <- runif(100,-3,3)
Z <- bsplines(x)
K <- attr(Z, "K")
attr(Z, "class")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
