#' plots the results of a
#' multinomial restricted unfolding
#'
#' @param x an object of type mru
#' @param dims which dimensions to visualize
#' @param ynames names of the response variables
#' @param ycol colour for representation of response variables
#' @param xcol colour for representation of predictor variables
#' @param ocol colour for representation of row objects
#' @param \dots additional arguments to be passed.
#' @return Plot of the results obtained from mru
#'
#' @examples
#' data(dataExample_mru)
#' y = as.matrix(dataExample_mru[1:20 , 1])
#' X = as.matrix(dataExample_mru[1:20 , 2:6])
#' output = mru(y = y, X = X, S = 2)
#' plot(output)
#'
#' @import ggforce
#' @import ggplot2
#' @import ggrepel
#'
#' @export
plot.mru = function(x, dims = c(1,2), ynames = NULL, ycol = "darkgreen", xcol = "lightskyblue", ocol = "grey",...)
{

  object = x

  # retrieving infomation from object
  U = object$U[, dims]
  V = object$V[, dims]
  X = object$X

  if(is.null(ynames)) ynames = object$ynames
  xnames = object$xnames

  UU = as.data.frame(U)
  VV = as.data.frame(V)
  colnames(UU) = colnames(VV) = c("dim1", "dim2")


  P = ncol(X)
  B = object$B[ , dims]
  Xo = object$Xoriginal

  # for solid line
  MCx1 <- data.frame(labs=character(),
                     varx = integer(),
                     dim1 = double(),
                     dim2 = double(), stringsAsFactors=FALSE)
  # for markers
  MCx2 <- data.frame(labs=character(),
                     varx = integer(),
                     dim1 = double(),
                     dim2 = double(), stringsAsFactors=FALSE)

  ll = 0
  lll = 0
  for(p in 1:P){
    b = matrix(B[p , ], 2, 1)
    # solid line
    minx = min(Xo[, p])
    maxx = max(Xo[, p])
    m.x1 = c(minx,maxx)
    markers1 = matrix((m.x1 - object$mx[p])/object$sdx[p], 2, 1)
    markerscoord1 = outer(markers1, b) # markers1 %*% t(b %*% solve(t(b) %*% b))
    MCx1[(ll + 1): (ll + 2), 1] = paste0(c("min", "max"), p)
    MCx1[(ll + 1): (ll + 2), 2] = p
    MCx1[(ll + 1): (ll + 2), 3:4] = markerscoord1
    ll = ll + 2
    # markers
    m.x2 = pretty(Xo[, p])
    m.x2 = m.x2[which(m.x2 > minx & m.x2 < maxx)]
    l.m = length(m.x2)
    markers2 = matrix((m.x2 - object$mx[p])/object$sdx[p], l.m, 1)
    markerscoord2 = outer(markers2, b) # markers2 %*% t(b %*% solve(t(b) %*% b))
    MCx2[(lll + 1): (lll + l.m), 1] = paste(m.x2)
    MCx2[(lll + 1): (lll + l.m), 2] = p
    MCx2[(lll + 1): (lll + l.m), 3:4] = markerscoord2
    lll = lll + l.m
  } # loop p


  plt = ggplot() +
    geom_point(data = VV, aes_string(x = 'dim1', y = 'dim2'), colour = ycol, size = 3) +
    geom_point(data = UU, aes_string(x = 'dim1', y = 'dim2'), col = ocol) +
    geom_text_repel(data = VV, aes_string(x = 'dim1', y = 'dim2', label = 'ynames', family = '"mono"')) +
    xlab("Dimension 1") +
    ylab("Dimension 2")

  ######################################################
  # variable axes with ticks and markers for predictors
  ######################################################
  plt = plt + geom_abline(intercept = 0, slope = B[,2]/B[,1], colour = xcol, linetype = 3) +
    geom_line(data = MCx1, aes_string(x = 'dim1', y = 'dim2', group = 'varx'), col = xcol, size = 1) +
    geom_point(data = MCx2, aes_string(x = 'dim1', y = 'dim2'), col = xcol) +
    geom_text(data = MCx2, aes_string(x = 'dim1', y = 'dim2', label = 'labs'), nudge_y = -0.08, size = 1.5)

  ######################################################
  # variable labels
  ######################################################
  a = ceiling(max(abs(c(ggplot_build(plt)$layout$panel_scales_x[[1]]$range$range, ggplot_build(plt)$layout$panel_scales_y[[1]]$range$range))))

  idx = apply(abs(B), 1, which.max)
  t = s = rep(NA,P)
  for(pp in 1:P){
    t[pp] = (a * 1.1)/(abs(B[pp,idx[pp]])) * B[pp,-idx[pp]]
    s[pp] = sign(B[pp,idx[pp]])
  }
  CC = cbind(idx, t, s)
  rownames(CC) = xnames

  bottom = which(CC[, "idx"] == 2 & CC[, "s"] == -1)
  top =  which(CC[, "idx"] == 2 & CC[, "s"] == 1)
  right = which(CC[, "idx"] == 1 & CC[, "s"] == 1)
  left = which(CC[, "idx"] == 1 & CC[, "s"] == -1)

  if(length(CC[top, "t"])==0){
    plt = plt + scale_x_continuous(limits = c(-a,a), breaks = CC[bottom, "t"], labels = rownames(CC)[bottom],
                                   sec.axis = sec_axis(trans ~ ., breaks = 0, labels = ""))
  }else{
    plt = plt + scale_x_continuous(limits = c(-a,a), breaks = CC[bottom, "t"], labels = rownames(CC)[bottom],
                                   sec.axis = sec_axis(trans ~ ., breaks = CC[top, "t"], labels = rownames(CC)[top]))
  }


  if(length(CC[right, "t"])==0){
    plt = plt + scale_y_continuous(limits = c(-a,a), breaks = CC[left, "t"], labels = rownames(CC)[left],
                                   sec.axis = sec_axis(trans ~ ., breaks = 0, labels = ""))
  }else{
    plt = plt + scale_y_continuous(limits = c(-a,a), breaks = CC[left, "t"], labels = rownames(CC)[left],
                                   sec.axis = sec_axis(trans ~ ., breaks = CC[right, "t"], labels = rownames(CC)[right]))
  }

  # plt = plt + scale_x_continuous(limits = c(-a,a), breaks = CC[bottom, "t"], labels = rownames(CC)[bottom],
  #                                sec.axis = sec_axis(trans ~ ., breaks = CC[top, "t"], labels = rownames(CC)[top]))
  # plt = plt + scale_y_continuous(limits = c(-a,a), breaks = CC[left, "t"], labels = rownames(CC)[left],
  #                                sec.axis = sec_axis(trans ~ ., breaks = CC[right, "t"], labels = rownames(CC)[right]))

  plt = plt + theme(axis.line = element_line(colour = "black"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank())

  plt = plt + coord_fixed()
  suppressWarnings(print(plt))

}
