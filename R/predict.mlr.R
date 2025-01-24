#'  The function predict.mlr makes predictions for a test/validation set
#'  based on a fitted mlr model
#'
#' @param object An \code{mlr} object
#' @param newX An N by P matrix with predictor variables for a test/validation set
#' @param \dots additional arguments to be passed.
#'
#' @return This function returns an object of the class \code{p.mlr} with components:
#' \item{Ghat}{Predicted values (probabilities) for the test set}
#'
#' @examples
#' \dontrun{
#' data(dataExample_mru)
#' y = as.matrix(dataExample_mru[ , 1])
#' X = as.matrix(dataExample_mru[ , 2:6])
#' output = mlr(y = y, X = X, base = 1)
#' preds = predict(output, newX = X[1:4, ])
#' }
#'
#'
#' @export

predict.mlr = function(object, newX, ...){
  X = scale(newX, center = object$mx, scale = object$sdx)
  ones.n = matrix(1, nrow(X), 1)
  theta = ones.n %*% t(object$m) + X %*% object$A
  Ghat = exp(theta) / rowSums(exp(theta))
  return(Ghat)
}
