% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/procx.R
\name{procx}
\alias{procx}
\title{Helper function for pre-processing the predictors}
\usage{
procx(X)
}
\arguments{
\item{X}{An N by P matrix with predictor variables}
}
\value{
Xoriginal

dichotomous indicator which predictor variables are dichotomous

X standardized matrix

mx averages of original variables

sdx standard deviation of original variables
}
\description{
Continuous predictor variables are standardized - mean zero, standard deviation one.
Categorical predictor variables are represented as dummy (0, 1) variables.
}
