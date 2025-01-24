#' Plots a Cumulative Logistic MDU model
#'
#' @param x an object of type clmdu
#' @param dims which dimensions to visualize
#' @param circles which circles to visualize
#' @param ycol colour for representation of response variables
#' @param xcol colour for representation of predictor variables
#' @param ocol colour for representation of row objects
#' @param markersize size of points
#' @param labelsize size of labels
#' @param \dots additional arguments to be passed.
#' @return Plot of the results obtained from clmdu
#'
#' @examples
#' \dontrun{
#' data(dataExample_clmdu)
#' Y = as.matrix(dataExample_clmdu[ , 1:8])
#' X = as.matrix(dataExample_clmdu[ , 9:13])
#' # unsupervised
#' output = clmdu(Y = Y, S = 2)
#' plot(output)
#' }
#'
#'
#' @import ggforce
#' @import ggplot2
#' @import ggrepel
#'
#'
#' @export

plot.clmdu = function(x, dims = c(1,2), circles = seq(1,R),
                      ycol = "darkgreen", xcol = "darkblue", ocol = "grey",
                      markersize = 2.5, labelsize = 3, ...)
{

  object = x
  ccol = ycol

  # retrieving information from object
  m = object$m
  U = object$U[, dims]
  V = object$V[, dims]
  R = nrow(V)
  X = object$X
  isx = !is.null(X)

  ynames = object$ynames
  xnames = object$xnames

  UU = as.data.frame(U)
  VV = as.data.frame(V)
  colnames(UU) = colnames(VV) = c("dim1", "dim2")
  rownames(VV) = ynames

  if(!is.null(circles)){
    a = -m[[circles[1]]]
    # circles1 = cbind(circles[1], outer( rep(1, length( m[[circles[1]]][m[[circles[1]]] > 0] ) ) , V[circles[1], ]), m[[circles[1]]][m[[circles[1]]] > 0] )
    circles1 = cbind(circles[1], outer( rep(1, length( a[a > 0] ) ) , V[circles[1], ]), a[a > 0] )
    if(length(circles)>1){
      for(r in circles[-1]){
        a = -1 * m[[r]]
        circles1 = rbind(circles1,
                         cbind(r, outer( rep(1, length( a[a > 0] ) ) , V[r, ]), a[a > 0]  )
        )
      }
    }

    colnames(circles1) = c("var", "x0", "y0", "r")
    cnames = rownames(circles1)
    circles2 = data.frame(circles1)
    # for labels
    circles3 = circles2
    circles3$x0 = circles2$x0 - sign(circles2$x0) * circles2$r
  }

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
  # baseplot with points for objects and response variables
  ######################################################
  plt = ggplot() +
    geom_point(data = UU, aes(x = .data$dim1, y = .data$dim2), col = ocol, alpha = 0.1) +
    # geom_point(data = UU, aes(x = .data$dim1, y = .data$dim2), col = ocol) +
    geom_point(data = VV, aes(x = .data$dim1, y = .data$dim2), colour = ycol, size = 5) 

  # margins <- c("l" = ggplot_build(plt)$layout$panel_scales_x[[1]]$range$range[1] - .1,
  #              "r" = ggplot_build(plt)$layout$panel_scales_x[[1]]$range$range[2] + .1,
  #              "b" = ggplot_build(plt)$layout$panel_scales_y[[1]]$range$range[1] - .1,
  #              "t" = ggplot_build(plt)$layout$panel_scales_y[[1]]$range$range[2] + .1)

  ######################################################
  # add circles representin thresholds
  ######################################################
  if(!is.null(circles)){
    plt = plt + geom_circle(data = circles2,
                          aes( x0 = .data$x0, y0 = .data$y0, r = r),
                          colour = ccol, fill= ccol,
                          alpha = 0.05, linetype = "dotted", show.legend = FALSE) +
    geom_text(data = circles3, aes(x = .data$x0, y = .data$y0, label = cnames), nudge_x = 0.1, size = 2)
  }

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
                                             family = 'mono', fontface = 'bold'), size = 5)
  }

  ######################################################
  # response variable labels
  ######################################################
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
      geom_label_repel(data = df2,
                       aes(x = .data$dim1, y = .data$dim2, label = .data$var, colour = .data$type, family = "mono", fontface = "bold"),
                       size = labelsize, show.legend = FALSE) +
      scale_color_manual(values = c("0" = xcol, "1" = ycol)) 
  }
  
  plt = plt + coord_fixed() + theme_lmda()

  suppressWarnings(print(plt))

  return(plt)
}
