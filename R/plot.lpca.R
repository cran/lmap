#' Plots a Logistic PCA Model
#'
#' @param x an object of type lpca
#' @param dims which dimensions to visualize
#' @param type either H (hybrid), I (inner product/pca), or D (distance/melodic)
#' @param pmarkers vector or list of length R (the number of response variables) with values between 0 and 1 for markers on the response variable axes.
#' @param ycol colour for representation of response variables
#' @param xcol colour for representation of predictor variables
#' @param ocol colour for representation of row objects
#' @param \dots additional arguments to be passed.
#' @return Plot of the results obtained from lpca
#'
#' @examples
#' \dontrun{
#' data(dataExample_lpca)
#' Y = as.matrix(dataExample_lpca[, 1:8])
#' X = as.matrix(dataExample_lpca[, 9:13])
#' # unsupervised
#' output = lpca(Y = Y, S = 2)
#' plot(output)
#' }
#'
#' @import ggforce
#' @import ggplot2
#' @import ggrepel
#'
#' @export
plot.lpca <- function(x, dims = c(1,2), type = "H",
                      pmarkers = seq(0.1, 0.9, by = 0.1),
                      ycol = "darkgreen", xcol = "darkblue", ocol = "grey",...)
{

  object = x

  if(type == "H"){
    plt = plot.lpcah(object, dims = dims, pmarkers = pmarkers, ycol = ycol, xcol = xcol, ocol = ocol)
  }
  if(type == "I"){
    plt = plot.lpca1(object, dims = dims, pmarkers = pmarkers, ycol = ycol, xcol = xcol, ocol = ocol)
  }
  if(type == "D"){
    plt = plot.lpca2(object, dims = dims, ycol = ycol, xcol = xcol, ocol = ocol)
  }
  suppressWarnings(print(plt))

  return(plt)
}





