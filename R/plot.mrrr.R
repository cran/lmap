#' Plots a Multinomial Reduced Rank Model
#'
#'
#' @param x an object of type mrrr
#' @param dims which dimensions to visualize
#' @param ycol colour for representation of response variables
#' @param xcol colour for representation of predictor variables
#' @param ocol colour for representation of row objects
#' @param markersize determines the size of the markera
#' @param labelsize determines the size of the labels
#' @param \dots additional arguments to be passed.
#' @return Plot of the results obtained from mrrr
#'
#' @examples
#' \dontrun{
#' data(dataExample_mru)
#' y = as.matrix(dataExample_mru[ , 1])
#' X = as.matrix(dataExample_mru[ , 2:6])
#' output = mru(y = y, X = X, S = 2)
#' plot(output)
#' }
#'
#' @import ggforce
#' @import ggplot2
#' @import ggrepel
#'
#' @export

plot.mrrr = function(x, dims = c(1,2),
                     ycol = "darkgreen", xcol = "darkblue", ocol = "grey",
                     markersize = 2.5, labelsize = 3, ...){
  object = x

  # library(ggplot2)
  # library(ggrepel)
  # source("~/surfdrive/LogitMDA/lmap-package/new/R/theme_lmda.R")
  # source("~/surfdrive/LogitMDA/lmap-package/new/R/make.df.for.varlabels.R")
  # source("~/surfdrive/LogitMDA/lmap-package/new/R/make.dfs.for.X.R")

  # retrieving information from object
  m = object$m
  U = object$U[, dims]
  V = object$V[, dims]
  X = object$X

  ynames = object$ynames
  xnames = object$xnames

  UU = as.data.frame(U)
  VV = as.data.frame(V)
  colnames(UU) = colnames(VV) = c("dim1", "dim2")

  P = ncol(X)
  C = nrow(V)
  B = object$B[ , dims]
  Xo = object$Xoriginal

  ######################################################
  # retrieve information for predictor variables variable axes
  ######################################################
  dfxs = make.dfs.for.X(Xo, P, B, object$xnames, object$mx, object$sdx)
  MCx1 = dfxs$MCx1
  MCx2 = dfxs$MCx2
  MCx3 = dfxs$MCx3
  dichotomous = dfxs$dichotomous

  ######################################################
  # retrieve information for response variable axes
  ######################################################

  VV1 = outer(object$V[, dims[1]], object$V[ , dims[1]], "-")
  VV2 = outer(object$V[, dims[2]], object$V[ , dims[2]], "-")
  diflabs = vector()
  difm = vector()
  lolist = list()
  t = 0
  for(i in 1:(C-1)){
    for(j in (i+1):C){
      diflabs = c(diflabs, paste0(ynames[i], "-", ynames[j]))
      difm = c(difm, object$m[i] - object$m[j])
      t = t + 1
      lolist[[t]] = pretty(log( object$Ghat[, i] / object$Ghat[, j] ) )
    }
  }

  VV = cbind(-VV1[lower.tri(VV1)], -VV2[lower.tri(VV2)])
  CC = nrow(VV)

  MCy <- data.frame(labs=character(),
                    vary = integer(),
                    dim1 = double(),
                    dim2 = double(), stringsAsFactors=FALSE)
  ll = 0
  for(c in 1:CC){
    lo = lolist[[c]]
    l.m = length(lo)
    markerlabs = paste(lo)
    markers = matrix(lo - difm[c], length(lo), 1)
    v = matrix(VV[c, dims], nrow = 2, ncol = 1)
    markerscoord = markers %*% t(v %*% solve(t(v) %*% v))
    MCy[(ll + 1): (ll + l.m), 1] = markerlabs
    MCy[(ll + 1): (ll + l.m), 2] = c
    MCy[(ll + 1): (ll + l.m), 3:4] = markerscoord
    ll = ll + l.m
  }

  ######################################################
  # plotting - objects
  ######################################################

  plt = ggplot() +
    geom_point(data = UU, aes(x = .data$dim1, y = .data$dim2), col = ocol) +
    geom_abline(intercept = 0, slope = VV[, 2]/VV[, 1], colour = ycol) +
    geom_abline(intercept = 0, slope = B[!dichotomous, 2]/B[!dichotomous, 1],
                colour = xcol, linetype = "dotted") +
    geom_line(data = MCx1, aes(x = .data$dim1, y = .data$dim2, group =.data$varx), col = xcol) +
    geom_label(data = MCx2, aes(x = .data$dim1, y = .data$dim2, label = .data$labs),
               fill = xcol, fontface = "bold", color = "white", size = markersize) +
    geom_point(data = MCx3, aes(x = .data$dim1, y = .data$dim2), col = xcol, size = 4) +
    geom_text(data = MCx3[-1, ], aes(x = .data$dim1, y = .data$dim2, label = labs,
                                     family = 'mono', fontface = 'bold'), nudge_y = 0.2, size = 3) +
    geom_label(data = MCy, aes(x = .data$dim1, y = .data$dim2, label = labs),
               fill = ycol, fontface = "bold", color = "white", size = markersize)

  margins <- c("l" = ggplot_build(plt)$layout$panel_scales_x[[1]]$range$range[1],
               "r" = ggplot_build(plt)$layout$panel_scales_x[[1]]$range$range[2],
               "b" = ggplot_build(plt)$layout$panel_scales_y[[1]]$range$range[1],
               "t" = ggplot_build(plt)$layout$panel_scales_y[[1]]$range$range[2])


  ######################################################
  # variable labels
  ######################################################

  PP = nrow(B[!dichotomous, ])
  df2 = make.df.for.varlabels(BV = rbind(B[!dichotomous, ], VV),
                              names = c(xnames[!dichotomous], diflabs),
                              margins = margins, P = PP, R = CC)
  df2$type = factor(df2$type)

  plt = plt +
    geom_label(data = df2,
                     aes(x = .data$dim1, y = .data$dim2, label = .data$var, colour = .data$type, family = "mono", fontface = "bold"),
                     size = labelsize, show.legend = FALSE) +
    scale_color_manual(values = c("0" = xcol, "1" = ycol)) +
    coord_fixed() +
    theme_lmda()

  suppressWarnings(print(plt))

  return(plt)
}
