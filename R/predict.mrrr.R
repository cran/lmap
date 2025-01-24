#'  The function predict.mrrr makes predictions for a test/validation set
#'  based on a fitted mrrr model
#'
#' @param object An \code{mrrr} object
#' @param newX An N by P matrix with predictor variables for a test/validation set
#' @param \dots additional arguments to be passed.
#'
#' @return This function returns an object of the class \code{p.mru} with components:
#' \item{Ghat}{Predicted values for the test set}
#'
#' @examples
#' \dontrun{
#' data(dataExample_lpca)
#' Y = as.matrix(dataExample_mru[-c(1:20) , 1:8])
#' X = as.matrix(dataExample_mru[-c(1:20) , 9:13])
#' newY = as.matrix(dataExample_mru[1:20 , 1:8])
#' newX = as.matrix(dataExample_mru[1:20 , 9:13])
#' output = mrrr(Y = Y, X = X, S = 2)
#' preds = predict(output, newX = newX)
#' }
#'
#'
#' @export

predict.mrrr = function(object, newX,...){
  X = scale(newX, center = object$mx, scale = object$sdx)
  ones.n = matrix(1, nrow(X), 1)
  theta = ones.n %*% t(object$m) + X %*% object$B %*% t(object$V)
  Ghat = exp(theta) / rowSums(exp(theta))
  return(Ghat)
}
