#' Summarizing Multinomial Restricted Unfolding Model
#' 											   
#' The function summary.mru gives a summary from an object from mru()
#'
#' @param object An object resulting from mru
#' @param \dots additional arguments to be passed.
#'
#' @return Summary of the results obtained from mru
#'
#' @export
summary.mru = function(object,...){

  cat("\n")
  cat("Call:", "\n")
  print(object$call)
  cat("\n")
  cat("Fitted multinomial restricted unfolding model", "\n")
  cat("Residual deviance:", object$deviance, "\n")
  cat("Number of fitted parameters:", object$npar, "\n" )
  cat("AIC:", object$AIC, "\n" )
  cat("BIC:", object$BIC, "\n" )
  cat("\n")
  cat("Bx coefficients:", "\n")
  print(round(object$Bx, digits = 2))
  cat("\n")
  if(!is.null(object$Bz)){
    cat("Bz coefficients:", "\n")
    print(round(object$Bz, digits = 2))
  }
  cat("V coefficients:", "\n")
  print(round(object$V, digits = 2))
  cat("\n")
}
