#' Plots a Cumulative Logistic PCA model
#'
#' @param x an object of type clpca
#' @param dims which dimensions to visualize
#' @param ycol colour for representation of response variables
#' @param xcol colour for representation of predictor variables
#' @param ocol colour for representation of row objects
#' @param markersize size of points
#' @param labelsize size of labels
#' @param \dots additional arguments to be passed.
#'
#' @return Plot of the results obtained from clpca
#' @examples
#' \dontrun{
#' data(dataExample_clpca)
#' Y<-as.matrix(dataExample_clpca[,5:8])
#' X<-as.matrix(dataExample_clpca[,1:4])
#' out = clpca(Y, X)
#' plot(out)
#' }
#'
#' @import ggplot2
#' @export
plot.clpca <- function(x, dims = c(1,2),
                       ycol = "darkgreen", xcol = "darkblue", ocol = "grey",
                       markersize = 2.5, labelsize = 3, ...)
{

  object<-x
  ######################################################
  # retrieve information from object
  ######################################################
  Y = object$Y

  U = as.data.frame(object$U[ , dims])
  N = nrow(U)
  colnames(U) = c("dim1", "dim2")

  V = object$V[ , dims]
  VV = as.data.frame(V)
  R = nrow(V)
  colnames(VV) = c("dim1", "dim2")
  rownames(VV) = object$ynames

  ######################################################
  # retrieve information for response variables variable axes
  ######################################################

  MCy <- data.frame(labs=character(),
                    vary = integer(),
                    dim1 = double(),
                    dim2 = double(), stringsAsFactors=FALSE)
  ll = 0

  for(r in 1:R){
    m = object$m[[r]]
    l.m = length(m)
    markers = matrix(m, ncol = 1)
    v = matrix(V[r, ], nrow = 2, ncol = 1)
    markerscoord = markers %*% t(v %*% solve(t(v) %*% v))
    markerlabs = names(m)
    MCy[(ll + 1): (ll + l.m), 1] = markerlabs
    MCy[(ll + 1): (ll + l.m), 2] = r
    MCy[(ll + 1): (ll + l.m), 3:4] = markerscoord
    ll = ll + l.m
  }

  ######################################################
  # retrieve information for predictor variables variable axes
  ######################################################
  X = object$X
  isx = !is.null(X)
  if(isx){
    P = ncol(X)
    B = object$B[ , dims]
    Xo = object$Xoriginal

    dfxs = make.dfs.for.X(Xo, P, B, object$xnames, object$mx, object$sdx)
    MCx1 = dfxs$MCx1
    MCx2 = dfxs$MCx2
    MCx3 = dfxs$MCx3
    dichotomous = dfxs$dichotomous
  } #isx

  ######################################################
  # plotting - objects
  ######################################################
  plt = ggplot() +
    geom_point(data = U, aes(x = .data$dim1, y = .data$dim2), colour = ocol) +
    xlab(paste("Dimension", dims[1])) +
    ylab(paste("Dimension", dims[2]))

  margins <- c("l" = ggplot_build(plt)$layout$panel_scales_x[[1]]$range$range[1] - .1,
               "r" = ggplot_build(plt)$layout$panel_scales_x[[1]]$range$range[2] + .1,
               "b" = ggplot_build(plt)$layout$panel_scales_y[[1]]$range$range[1] - .1,
               "t" = ggplot_build(plt)$layout$panel_scales_y[[1]]$range$range[2] + .1)

  ######################################################
  # variable axes with ticks and markers for predictors
  ######################################################
  if(isx){
    plt = plt +
      geom_abline(intercept = 0, slope = B[!dichotomous, 2]/B[!dichotomous,1], colour = xcol, linetype = 3) +
      geom_line(data = MCx1, aes(x = .data$dim1, y = .data$dim2, group = .data$varx), col = xcol, linewidth = 1) +
      geom_label(data = MCx2, aes(x = .data$dim1, y = .data$dim2, label = .data$labs),
                 fill = xcol, fontface = "bold", color = "white", size = markersize) +
      geom_point(data = MCx3, aes(x = .data$dim1, y = .data$dim2), col = xcol, size = 4) +
      geom_text_repel(data = MCx3[-1, ], aes(x = .data$dim1, y = .data$dim2, label = labs,
                                             family = 'mono', fontface = 'bold'), size = 3)
  }

  ######################################################
  # variable axes with ticks and markers for responses
  ######################################################
  plt = plt +
    geom_abline(intercept = 0, slope = V[,2]/V[,1], colour = ycol) +
    geom_label(data = MCy, aes(x = .data$dim1, y = .data$dim2, label = labs),
               fill = ycol, fontface = "bold", color = "white", size = markersize)

  ######################################################
  # variable labels
  ######################################################
  margins <- c("l" = ggplot_build(plt)$layout$panel_scales_x[[1]]$range$range[1],
               "r" = ggplot_build(plt)$layout$panel_scales_x[[1]]$range$range[2],
               "b" = ggplot_build(plt)$layout$panel_scales_y[[1]]$range$range[1],
               "t" = ggplot_build(plt)$layout$panel_scales_y[[1]]$range$range[2])

  if(isx){
    BV = rbind(B[!dichotomous, ], V)
    PP = nrow(B[!dichotomous, ])
    CC = nrow(V)
    names = c(object$xnames[!dichotomous], object$ynames)
  }
  else{
    BV = V
    PP = 0
    CC = nrow(V)
    names = object$ynames
  }

  df2 = make.df.for.varlabels(BV = BV, names = names,
                              margins = margins, P = PP, R = CC)
  df2$type = factor(df2$type)

  plt = plt +
    geom_label_repel(data = df2,
                     aes(x = .data$dim1, y = .data$dim2, label = .data$var, colour = .data$type, family = "mono", fontface = "bold"),
                     size = labelsize, show.legend = FALSE) +
    scale_color_manual(values = c("0" = xcol, "1" = ycol)) +
    coord_fixed() +
    theme_lmda()

  suppressWarnings(print(plt))
  return(plt)
}
