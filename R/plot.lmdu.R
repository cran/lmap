#' Plots a Logistic MDU model
#'
#' @param x an object of type lmdu
#' @param dims which dimensions to visualize
#' @param ycol colour for representation of response variables
#' @param xcol colour for representation of predictor variables
#' @param ocol colour for representation of row objects
#' @param markersize size of points
#' @param labelsize size of labels
#' @param \dots additional arguments to be passed.
#' @return Plot of the results obtained from lmdu
#'
#' @examples
#' \dontrun{
#' data(dataExample_lmdu)
#' Y = as.matrix(dataExample_lmdu[ , 1:8])
#' X = as.matrix(dataExample_lmdu[ , 9:13])
#' # unsupervised
#' output = lmdu(Y = Y, S = 2)
#' plot(output)
#' }
#'
#' @import ggforce
#' @import ggplot2
#' @import ggrepel
#'
#' @export
plot.lmdu = function(x, dims = c(1,2),
                     ycol = "darkgreen", xcol = "darkblue", ocol = "grey",
                     markersize = 2.5, labelsize = 3, ...)
{

  object = x

  # retrieving information from object
  m = object$m
  U = object$U[, dims]
  V = object$V[, dims]
  X = object$X
  isx = !is.null(X)

  ynames = object$ynames
  xnames = object$xnames

  UU = as.data.frame(U)
  VV = as.data.frame(V)
  colnames(UU) = colnames(VV) = c("dim1", "dim2")

  radius = rep(NA,nrow(V))
  for(r in 1:nrow(V)){radius[r] = ifelse(m[r]>0, m[r] , 0)}
  circles = VV
  circles[,3] = radius
  colnames(circles) = c("x0","y0","r")

  if(isx){
    P = ncol(X)
    B = object$B[ , dims]
    Xo = object$Xoriginal

    dfxs = make.dfs.for.X(Xo, P, B, xnames, object$mx, object$sdx)
    MCx1 = dfxs$MCx1
    MCx2 = dfxs$MCx2
    MCx3 = dfxs$MCx3
    dichotomous = dfxs$dichotomous

  } #isx


  plt = ggplot() +
    geom_point(data = UU, aes(x = .data$dim1, y = .data$dim2), col = ocol) +
    geom_point(data = VV, aes(x = .data$dim1, y = .data$dim2), colour = ycol, size = 5) +
    geom_circle(data=circles, aes(x0 = .data$x0, y0 = .data$y0, r = r), colour = ycol, fill= ycol, alpha = 0.05, show.legend = FALSE) +
    xlab(paste("Dimension", dims[1])) +
    ylab(paste("Dimension", dims[2]))

  # margins <- c("l" = ggplot_build(plt)$layout$panel_scales_x[[1]]$range$range[1] - .1,
  #              "r" = ggplot_build(plt)$layout$panel_scales_x[[1]]$range$range[2] + .1,
  #              "b" = ggplot_build(plt)$layout$panel_scales_y[[1]]$range$range[1] - .1,
  #              "t" = ggplot_build(plt)$layout$panel_scales_y[[1]]$range$range[2] + .1)

  ######################################################
  # variable axes with ticks and markers for predictors
  ######################################################
  if(isx){
    plt = plt +
      geom_abline(intercept = 0, slope = B[!dichotomous, 2]/B[!dichotomous,1], colour = xcol, linetype = 3) +
      geom_line(data = MCx1, aes(x = .data$dim1, y = .data$dim2, group = .data$varx), col = xcol) +
      geom_label(data = MCx2, aes(x = .data$dim1, y = .data$dim2, label = .data$labs),
                 fill = xcol, fontface = "bold", color = "white", size = markersize) +
      geom_point(data = MCx3, aes(x = .data$dim1, y = .data$dim2), col = xcol, size = 4) +
      geom_text_repel(data = MCx3[-1, ], aes(x = .data$dim1, y = .data$dim2, label = labs,
                                             family = 'mono', fontface = 'bold'), size = 3)
  }

  rownamesVV<-rownames(VV)
  plt = plt + geom_text_repel(data = VV, aes(x = .data$dim1, y = .data$dim2, label = rownamesVV,
                                             family = 'mono', fontface = 'bold'), col = "black", size = 5)
  ######################################################
  # variable labels
  ######################################################

  if(isx){
    margins <- c("l" = ggplot_build(plt)$layout$panel_scales_x[[1]]$range$range[1],
                 "r" = ggplot_build(plt)$layout$panel_scales_x[[1]]$range$range[2],
                 "b" = ggplot_build(plt)$layout$panel_scales_y[[1]]$range$range[1],
                 "t" = ggplot_build(plt)$layout$panel_scales_y[[1]]$range$range[2])

    B = B[!dichotomous, ]
    PP = nrow(B)
    df2 = make.df.for.varlabels(BV = B, names = xnames[!dichotomous],
                                margins = margins, P = PP, R = 0)
    df2$type = factor(df2$type)

    plt = plt +
      geom_label(data = df2,
                       aes(x = .data$dim1, y = .data$dim2, label = .data$var, colour = .data$type, family = "mono", fontface = "bold"),
                       size = labelsize, show.legend = FALSE) +
      scale_color_manual(values = c("0" = xcol, "1" = ycol))
  }

  plt = plt + coord_fixed() + theme_lmda()
  
  suppressWarnings(print(plt))

  return(plt)
}
