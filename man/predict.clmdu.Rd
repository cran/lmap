% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.clmdu.R
\name{predict.clmdu}
\alias{predict.clmdu}
\title{The function predict.clmdu makes predictions for a test/validation set
 based on a fitted cl restricted multidimensional unfolding model (clmdu with X)}
\usage{
\method{predict}{clmdu}(object, newX, newY = NULL, ...)
}
\arguments{
\item{object}{An \code{clmdu} object}

\item{newX}{An N by P matrix with predictor variables for a test/validation set}

\item{newY}{An N by R matrix with response variables  for a test/validation set}

\item{\dots}{additional arguments to be passed.}
}
\value{
This function returns an object of the class \code{predclpca} with components:
\item{Yhat}{Predicted values for the test set}
\item{devr}{Estimated prediction deviance for separate responses}
\item{devtot}{Estimated prediction deviance for all responses}
}
\description{
The function predict.clmdu makes predictions for a test/validation set
 based on a fitted cl restricted multidimensional unfolding model (clmdu with X)
}
\examples{
\dontrun{
data(dataExample_clpca)
Y = as.matrix(dataExample_clmdu[ , 1:8])
X = as.matrix(dataExample_clmdu[ , 9:13])
newY = as.matrix(dataExample_clmdu[1:20 , 1:8])
newX = as.matrix(dataExample_clmdu[1:20 , 9:13])
# supervised
output = clmdu(Y = Y, X = X, S = 2)
preds = predict(output, newX = newX, newY = newY)
}

}
