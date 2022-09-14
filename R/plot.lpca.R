#' plots the results of a
#' logistic principal component analysis (X = NULL)
#' logistic reduced rank regression (X != NULL)
#'
#' @param x an object of type lpca
#' @param dims which dimensions to visualize
#' @param type either pca or dist
#' @param ycol colour for representation of response variables
#' @param xcol colour for representation of predictor variables
#' @param ocol colour for representation of row objects
#' @param \dots additional arguments to be passed.
#' @return Plot of the results obtained from lpca
#'
#' @examples
#' data(dataExample_lpca)
#' Y = as.matrix(dataExample_lpca[1:20 , 1:8])
#' X = as.matrix(dataExample_lpca[1:20 , 9:13])
#' # unsupervised
#' output = lpca(Y = Y, S = 2)
#' plot(output)
#'
#'
#' @import ggforce
#' @import ggplot2
#' @import ggrepel
#'
#' @export
plot.lpca <- function(x, dims = c(1,2), type = "pca", ycol = "darkgreen", xcol = "lightskyblue", ocol = "grey",...)
{

  object = x

  if(type == "pca"){
    plt = plot.lpca1(object, dims = c(1,2), ycol = "darkgreen", xcol = "lightskyblue", ocol = "grey")
  }
  if(type == "dist"){
    plt = plot.lpca2(object, dims = c(1,2), ycol = "darkgreen", xcol = "lightskyblue", ocol = "grey")
  }
  suppressWarnings(print(plt))

  return(plt)
}





