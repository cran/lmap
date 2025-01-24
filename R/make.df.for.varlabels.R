#' Helper function for the plot functions
#'
#' Helper function for the plot functions to add variable labels 
#' to the predictor and response variable axes
#'
#' @param BV Concatention of matrices B (P x S) and V (R x S) or only V (R x S) from lpca
#' @param margins a vector of length four indicating the margins of the plot
#' @param names a vector with variable names
#' @param P an integer indicating the number of predictor variables
#' @param R an integer indicating the number of response variables
#' @return df with information for the placement of the variable labels
#'
#' @export
make.df.for.varlabels = function(BV, margins, names, P, R){
  margins = margins - sign(margins) * 0.1 # plaats labels net binnen marges
  beta <- BV[,2]/BV[,1] # slope van de variabele axes

  lab <- data.frame("xname" = names,
                    "b" = beta,
                    "Yleft" = beta * margins["l"],
                    "Yright" = beta * margins["r"])

  orientation = sign(BV[,1]) #sign of dim1 defines direction l-r
  lab$side =  c("left","right")[ as.numeric( BV[ , 1] > 0 ) + 1]
  lab$side[lab$Yleft < margins["b"] & orientation < 0 ] = "bottom"
  lab$side[lab$Yleft > margins["t"] & orientation < 0 ] = "top"
  lab$side[lab$Yright < margins["b"] & orientation > 0] = "bottom"
  lab$side[lab$Yright > margins["t"] & orientation > 0] = "top"

  lab$Y <- lab$X <- NA
  lab$X[lab$side == "bottom"] = (margins["b"] / beta[lab$side == "bottom"])
  lab$Y[lab$side == "bottom"] = margins["b"]

  lab$X[lab$side == "top"] = (margins["t"] / beta[lab$side == "top"])
  lab$Y[lab$side == "top"] = margins["t"]

  lab$X[lab$side == "left"] = margins["l"]
  lab$Y[lab$side == "left"] = margins["l"] * beta[lab$side == "left"]

  lab$X[lab$side == "right"] = margins["r"]
  lab$Y[lab$side == "right"] = margins["r"] * beta[lab$side == "right"]

  df2 = lab[, c("xname", "X", "Y")]
  df2 = cbind(df2, c(rep(0, P), rep(1, R)))
  colnames(df2) = c("var", "dim1", "dim2", "type")
  return(df2)
}
