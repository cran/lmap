#'  Plotting function for object of class trioscale
#'
#' @param x An object of class trioscale.
#' @param ycol colour for representation of response variables
#' @param xcol colour for representation of predictor variables
#' @param ocol colour for representation of row objects
#' @param markersize size of points
#' @param labelsize size of labels
#' @param classlabels List with plotting options for the labels of the Anchor points
#' @param s1 scaling factor for distance between points and log-ratio axes
#' @param s2 scaling factor for positioning class labels
#' @param s3 scaling factor for positioning variable lables
#' @param \dots additional arguments to be passed.
#'
#' @return This function returns an plot
#'
#' @examples
#' \dontrun{
#' out = trioscale(data)
#' plot.trioscale(out)
#' }
#'
#'
#' @export

plot.trioscale = function(x,
                          ycol = "darkgreen", xcol = "darkblue", ocol = "grey",
                          markersize = 2.5, labelsize = 3, classlabels = NULL,
                          s1 = 2.5,  s2 = 1.05, s3 = 1.15, ...)
{

  object<-x
  # ----------------------------------------------------------------
  # ----------------------------------------------------------------
  # helper functions
  # ----------------------------------------------------------------
  # ----------------------------------------------------------------
  ts.lraxes = function(X, s1, s2, classlabels){
    # for log ratio axes
    LRaxis <- data.frame(labs=character(),
                         varx = integer(),
                         dim1 = double(),
                         dim2 = double(),
                         colx = character(), stringsAsFactors=FALSE)

    # for ticks
    LRticks <- data.frame(labs=character(),
                          varx = integer(),
                          dim1 = double(),
                          dim2 = double(), stringsAsFactors=FALSE)

    # for labels of tick marks
    LRlabels <- data.frame(labs=character(),
                           varx = integer(),
                           dim1 = double(),
                           dim2 = double(), stringsAsFactors=FALSE)

    # for labels of Classes
    Manchor = data.frame(labs=character(),
                         dim1 = double(),
                         dim2 = double(), stringsAsFactors=FALSE)

    t1 = 0 # tellertje
    t2 = 0 # tellertje
    t3 = 0 # tellertje

    for (k in 1:3) {

      # Select orientation
      phi = c(0, -4, 4)[k] * pi / 6
      cs = cos(phi)
      sn = sin(phi)

      # Compute rotated coordinates
      M = matrix(c(cs, sn, -sn, cs), 2, 2)
      G = X %*% M
      # lower = s1 * min(G[,2])
      lower = min(G[,2]) - s1
      # dx = diff(range(G[,2]))/10
      ur = cbind(range(G[, 1]), lower) %*% t(M)

      # Log-ratio axes
      LRaxis[(t1 + 1): (t1 + 2), 1] = paste0(c("min", "max"), k)
      LRaxis[(t1 + 1): (t1 + 2), 2] = k
      LRaxis[(t1 + 1): (t1 + 2), 3:4] = ur
      t1 = t1 + 2

      # Find tick marks in range of data
      tc = pretty(G[, 1])
      rg = range(G[, 1])
      tc = tc[rg[1] <= tc & tc <= rg[2]]
      nt = length(tc)
      ur = cbind(rep(tc, each = 2), rep(c(lower, lower - s1), nt)) %*% t(M)

      LRticks[(t2 + 1): (t2 + 2*nt), 1] = paste0(1:nt, "/" , k)
      LRticks[(t2 + 1): (t2 + 2*nt), 2] = rep(seq(from = ((k-1) * t3 + 1), to = ((k-1) * t3 + nt)), each = 2)
      LRticks[(t2 + 1): (t2 + 2*nt), 3:4] = ur
      t2 = t2 + 2 * nt

      # Find tick labels in range of data
      ur = cbind(tc, lower) %*% t(M)
      LRlabels[(t3 + 1): (t3 + nt), 1] = paste(abs(tc))
      LRlabels[(t3 + 1): (t3 + nt), 2] = seq(from = ((k-1) * nt + 1), to = (k * nt))
      LRlabels[(t3 + 1): (t3 + nt), 3:4] = ur
      t3 = t3 + nt
    }

    # LRaxis$colx = as.factor(LRaxis$colx)
    # LRticks$colx = as.factor(LRticks$colx)
    # LRlabels$colx = as.factor(LRlabels$colx)

    # Labels for the three classes
    Manchor[1:3 , 1] = classlabels
    # Manchor[1:3 , 2:3] = (LRaxis[c(1,2,3), 3:4] + LRaxis[c(4,5,6), 3:4])/2
    Manchor[1 , 2] = LRaxis[4, 3]
    Manchor[1 , 3] = LRaxis[1, 4]
    Manchor[2 , 2] = LRaxis[5, 3]
    Manchor[2 , 3] = LRaxis[2, 4]
    Manchor[3 , 2] = (LRaxis[3, 3] + LRaxis[6, 3])/2
    Manchor[3 , 3] = max(LRaxis[3,4], LRaxis[6,4])
    Manchor[, 2:3] = Manchor[ , 2:3] * s2 # scaling factor ipv 1.3

    # create output object.
    output = list(
      LRaxis = LRaxis,
      LRticks = LRticks,
      LRlabels = LRlabels,
      LRanchor = Manchor,
      dens = NULL
    )

    return(output)
  }

  # ----------------------------------------------------------------
  # ----------------------------------------------------------------
  # unpacking object
  # ----------------------------------------------------------------
  # ----------------------------------------------------------------

  mlr.output = object$mlr
  y = object$mlr$y
  X = object$mlr$Xoriginal
  U = object$U
  Udf = object$Udf

  dfxs = object$dfxs
  dichotomous = dfxs$dichotomous
  MCx1 = dfxs$MCx1
  MCx1a = MCx1[seq(2, nrow(MCx1), by = 2), ]
  MCx1a[ , 3:4] = s3 * MCx1a[ , 3:4]
  MCx1a$labs = object$mlr$xnames[!dichotomous]
  MCx2 = dfxs$MCx2
  MCx3 = dfxs$MCx3

  if(is.null(classlabels)){
    classlabels = mlr.output$ynames
  }

  # ----------------------------------------------------------------
  # ----------------------------------------------------------------
  # SOME COMPUTATIONS for LOG-RATIO AXES and TICKS
  # ----------------------------------------------------------------
  # ----------------------------------------------------------------

  # computations for Log-Ratio axes
  output.lraxes = ts.lraxes(U, s1 = s1, s2 = s2, classlabels)
  LRaxis = output.lraxes$LRaxis
  LRlabels = output.lraxes$LRlabels
  LRanchor = output.lraxes$LRanchor
  # ----------------------------------------------------------------
  # ----------------------------------------------------------------
  # PLOTTING
  # ----------------------------------------------------------------
  # ----------------------------------------------------------------
  myplot = ggplot()

  # add points for samples
  myplot = myplot +
    geom_point(data = Udf, aes(x = .data$dim1, y = .data$dim2), col = ocol) +
    coord_fixed()

  # log-ratio axes
  myplot = myplot +
    geom_line(data = LRaxis, aes(x = .data$dim1, y = .data$dim2, group = .data$varx),
              col = ycol, linewidth = 1) +
    # geom_line(data = output.lraxes$LRticks, aes(x = dim1, y = dim2, group = varx), col = ycolor) +
    geom_label(data = LRlabels, aes(x = .data$dim1, y = .data$dim2, label = .data$labs),
               fill = ycol, fontface = "bold", color = "white", size = markersize) +
    geom_label(data = LRanchor, aes(x = .data$dim1, y = .data$dim2, label = .data$labs),
               size = labelsize, family = 'mono')

  # predictor variable axes
  myplot = myplot +
    geom_line(data = MCx1, aes(x = .data$dim1, y = .data$dim2, group = .data$varx),
              col = xcol, linewidth = 1) +
    geom_label(data = MCx2, aes(x = .data$dim1, y = .data$dim2, label = .data$labs),
               fill = xcol, fontface = "bold", color = "white", size = markersize) +
    geom_label(data = MCx1a, aes(x = .data$dim1, y = .data$dim2, label = .data$labs),
               size = labelsize, family = 'mono') +
    geom_point(data = MCx3, aes(x = .data$dim1, y = .data$dim2), col = xcol, size = 4) +
    geom_text_repel(data = MCx3[-1, ], aes(x = .data$dim1, y = .data$dim2, label = labs,
                                           family = 'mono', fontface = 'bold'), size = 3)

  # remove potential legends
  myplot = myplot + theme_lmda()

  # plot
  suppressWarnings(print(myplot))

  # ----------------------------------------------------------------
  # ----------------------------------------------------------------
  # OUTPUT OBJECT
  # ----------------------------------------------------------------
  # ----------------------------------------------------------------
  return(myplot)
}
