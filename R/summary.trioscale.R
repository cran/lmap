#' Summarizing TrioScale
#'
#' The function summary.trioscale gives a summary from the MLR in trioscale
#'
#' @param object An object resulting from trioscale
#' @param \dots additional arguments to be passed.
#'
#' @return Summary of the results obtained from trioscale
#'
#' @export

summary.trioscale = function(object,...){
  summary.mlr(object$mlr)
}
