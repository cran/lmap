#' Summarizing Multinomial Logistic Regression Model
#'
#' The function summary.mlr gives a summary from an object from mlr()
#'
#' @param object An object resulting from mlr
#' @param \dots additional arguments to be passed.
#'
#' @return Summary of the results obtained from mlr
#'
#' @export

summary.mlr = function(object,...){
  cat("\n")
  cat("Call:", "\n")
  print(object$call)
  cat("\n")
  cat("Fitted multinomial regression model", "\n")
  cat("Residual deviance:", object$deviance, "\n")
  cat("Number of fitted parameters:", object$npar, "\n" )
  cat("AIC:", object$AIC, "\n" )
  cat("BIC:", object$BIC, "\n" )
  cat("\n")
  cat("Intercepts:", "\n")
  print(round(object$m, digits = 2))
  cat("\n")
  cat("Regression coefficients:", "\n")
  print(round(object$A, digits = 2))
  cat("\n")
}
