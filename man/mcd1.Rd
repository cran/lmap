% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcd1.R
\name{mcd1}
\alias{mcd1}
\title{Multinomial Canonical Decomposition Model for Multivariate Binary Data}
\usage{
mcd1(
  X,
  Y,
  S = 2,
  Z = NULL,
  W = NULL,
  ord.z = 1,
  ord.m = R,
  trace = FALSE,
  maxiter = 65536,
  dcrit = 1e-06
)
}
\arguments{
\item{X}{An N by P matrix with predictor variables}

\item{Y}{An N times R binary matrix  .}

\item{S}{Positive number indicating the dimensionality of teh solution}

\item{Z}{design matrix for response}

\item{W}{design matrix for intercepts}

\item{ord.z}{if Z = NULL, the function creates Z having order ord.z}

\item{ord.m}{if W = NULL, the function creates W having order ord.m}

\item{trace}{whether progress information should be printed on the screen}

\item{maxiter}{maximum number of iterations}

\item{dcrit}{convergence criterion}
}
\value{
This function returns an object of the class \code{mcd} with components:
\item{call}{function call}
\item{Xoriginal}{Matrix X from input}
\item{X}{Scaled X matrix}
\item{mx}{Mean values of X}
\item{sdx}{Standard deviations of X}
\item{Y}{Matrix Y from input}
\item{pnames}{Variable names of profiles}
\item{xnames}{Variable names of predictors}
\item{znames}{Variable names of responses}
\item{Z}{Design matrix Z}
\item{W}{Design matrix W}
\item{G}{Profile indicator matrix G}
\item{m}{main effects}
\item{bm}{regression weights for main effects}
\item{Bx}{regression weights for X}
\item{Bz}{regression weights for Z}
\item{A}{regression weights (Bx Bz')}
\item{U}{matrix with coordinates for row-objects}
\item{V}{matrix with coordinates for column-objects}
\item{Ghat}{Estimated values of G}
\item{deviance}{value of the deviance at convergence}
\item{df}{number of paramters}
\item{AIC}{Akaike's informatoin criterion}
\item{iter}{number of main iterations from the MM algorithm}
\item{svd}{Singular value decomposition in last iteration}
}
\description{
The function mcd1 fits the multinomial canonical decomposition model to multivariate binary responses
i.e. a double constrained reduced rank multinomial logistic model
}
\examples{
\dontrun{
data(dataExample_lpca)
Y = as.matrix(dataExample_lpca[ , 1:5])
X = as.matrix(dataExample_lpca[ , 9:13])
#unsupervised
output = mcd1(X, Y, S = 2, ord.z = 2)
}

}
