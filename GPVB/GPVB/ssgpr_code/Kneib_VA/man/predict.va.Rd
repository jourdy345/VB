 \name{predict.va}
\alias{predict.va}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict method for Variational Approximations
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Predicts values based on a variational approximation object
}

\usage{
\method{predict}{va}(object, newdata=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.

  \arguments{
 \item{object}{
 object of va-type}
 \item{newdata}{
   data for which the prediction is to be done. The data has to be passed in form of matrix or dataframe. The colnames have to be the ones, to which the formula referes. If newdata=NULL the prediction of the va object is returned 
   (usually this is the prediction based on the estimation data set, if the pred-Argument was NULL in the estimation. 
    If there was a prediction-dataset, predict returns the prediction for this dataset)
 }
\item{...}{Additional arguments}
  }


\value{
The prediction is returned as a vector
}

\author{
Elisabeth Waldmann
}
\note{
Random effects are set to zero in the prediction.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
 