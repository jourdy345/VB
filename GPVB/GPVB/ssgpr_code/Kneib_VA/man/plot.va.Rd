 \name{plot.va}
\alias{plot.va}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots for VA
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Plots the different effects of a va object.
}

\usage{
\method{plot}{va}(x,map=NULL,ask=TRUE,linear=TRUE,bsplines=TRUE,spatial=TRUE,random=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.

  \arguments{
 \item{x}{
 object of va-type}
 \item{map}{
 if there is a spatial effect, the map has to be given as boundary file (i.e. .bnd)}
\item{ask}{
 if \code{ask=TRUE}, the user will be asked before the page changes to the following plot}
 \item{linear}{if \code{linear = FALSE} there will be no plot of the linear effects, even though there are linear effects in the model}
 \item{bsplines}{if \code{bsplines = FALSE}, there will be no plot of the nonlinear effects, even though there are nonlinear effects in the model}
 \item{spatial}{if \code{spatial = FALSE}, there will be no plot of the spatial effects, even though there are spatial effects in the model}
 \item{random}{if \code{random = FALSE}, there will be no plot of the random effects, even though there are random effects in the model}
 \item{...}{Additional arguments}
}




\details{
The function returns a plot of densities of the linear effects, the fitted curve of the centralized b-splines (plotted against the covariates), a map with colors refering to the spatial effects and barplots of the random effects. The map is plotted with the function \link{drawmap} of the package \pkg{BayesX}.
%%  ~~ If necessary, more details than the description above ~~
}

\author{
Elisabeth Waldmann \cr
\email{ewaldma@uni-goettingen.de}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
 