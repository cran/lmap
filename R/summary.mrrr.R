#' Summarizing Multinomial Reduced Rank Model
#'
#' The function summary.mrrr gives a summary from an object from mrrr()
#'
#' @param object An object resulting from mrrr
#' @param \dots additional arguments to be passed.
#'
#' @return Summary of the results obtained from mrrr
#'
#' @export

summary.mrrr = function(object,...){
  cat("\n")
  cat("Call:", "\n")
  print(object$call)
  cat("\n")
  cat("Fitted multinomial reduced rank model", "\n")
  cat("Residual deviance:", object$deviance, "\n")
  cat("Number of fitted parameters:", object$npar, "\n" )
  cat("AIC:", object$AIC, "\n" )
  cat("BIC:", object$BIC, "\n" )
  cat("\n")
  cat("Intercepts:", "\n")
  print(round(object$m, digits = 2))
  cat("\n")
  cat("B coefficients:", "\n")
  print(round(object$B, digits = 2))
  cat("\n")
  cat("V coefficients:", "\n")
  print(round(object$V, digits = 2))
  cat("\n")
  cat("Implied coefficients:", "\n")
  print(round(object$B %*% t(object$V), digits = 2))
  cat("\n")
}
