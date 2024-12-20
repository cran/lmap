% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.clpca.R
\name{predict.clpca}
\alias{predict.clpca}
\title{The function predict.clpca makes predictions for a test/validation set
 based on a fitted clrrr model (clpca with X)}
\usage{
\method{predict}{clpca}(object, newX, newY = NULL, ...)
}
\arguments{
\item{object}{An \code{clpca} object}

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
The function predict.clpca makes predictions for a test/validation set
 based on a fitted clrrr model (clpca with X)
}
\examples{
\dontrun{
data(dataExample_clpca)
Y = as.matrix(dataExample_clpca[ , 1:8])
X = as.matrix(dataExample_clpca[ , 9:13])
newY = as.matrix(dataExample_clpca[1:20 , 1:8])
newX = as.matrix(dataExample_clpca[1:20 , 9:13])
# supervised
output = clpca(Y = Y, X = X, S = 2)
preds = predict(output, newX = newX, newY = newY)
}

}
