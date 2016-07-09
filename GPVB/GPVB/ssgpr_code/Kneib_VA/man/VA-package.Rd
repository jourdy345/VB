\name{VA-package}
\alias{VA-package}
\alias{VA}
\docType{package}
\title{
Variational Approximations and its evaluations in geoadditive (quantile) regression

}
\description{
This package gives the possibility to estimate Bayesian additive (quantile) regression estimated via variational approximations.
}
\details{
\tabular{ll}{
Package: \tab VA\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-10-17\cr
License: 2.14.0
LazyLoad: \tab yes\cr
}
The main function of this package is \link{va}, a function which uses the theory of variational approximations to calculate the parameters for Bayesian geoadditive regression models for mean regression as well as for quantile regression. The other functions are or for data preperation for the different effects (\link{lin}, \link{bsplines}, \link{spatial})and (\link{risk})for model evaluation.
}
\references{
Yue, R.Y., Rue, H. \emph{Bayesian inference for additive mixed quantile regression models} Computational Statistics & Data Analysis, (2011),   55, 84-96.


Ormerod, J.T. and Wand, M.P.
\emph{Explaining Variational Approximations}.
The American Statistician, (2010), 64, 140-153. 

Wand, M. P.; Ormerod, J. T.; Padoan, S. A.; and Fruhwirth, R., \emph{Variational Bayes for Elaborate Distributions},
\url{http://ro.uow.edu.au/cssmwp/56} (working paper).


}
\author{
Elisabeth Waldmann \cr
University of Goettingen \cr
\email{ewaldma@uni-goettingen.de}\cr

with contributions from\cr
\cr
Thomas Kneib\cr
Universtity of Goettingen\cr
\url{http://www.uni-goettingen.de/de/264255.html}\cr

Udo Schroeder\cr
University of Oldenburg\cr
\url{http://www.uni-oldenburg.de/index/personen/?username=USchroeder1}\cr

Maintainer: Elisabeth Waldmann
}
\keyword{ Variational Approximation
Quantile Regression
Additive Regression}
\seealso{
\pkg{BayesX}
}
