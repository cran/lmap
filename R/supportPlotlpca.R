plot.lpca1 <- function(object, dims = c(1,2), pmarkers,
                       ycol = "darkgreen", xcol = "darkblue", ocol = "grey",
                       markersize = 2.5, labelsize = 3)
{
  # plots the results of a
  # logistic principal component analysis (X = NULL)
  # logistic reduced rank regression (X != NULL)
  #
  # @param object an object of type lpca
  # @param dims which dimensions to visualize
  # @param ycol colour for representation of response variables
  # @param xcol colour for representation of predictor variables
  # @param ocol colour for representation of row objects
  #
  
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
  
  if(!is.list(pmarkers)){
    pp = pmarkers
    pmarkers <- vector(mode = "list", length = R)
    for(r in 1:R){pmarkers[[r]] = pp}
  }
  
  for(r in 1:R){
    probs = pmarkers[[r]]
    lo = log(probs/(1 - probs))
    l.m = length(lo)
    markerlabs = paste(probs)
    
    markers = matrix(lo - object$m[r], length(lo), 1)
    v = matrix(V[r, dims], nrow = 2, ncol = 1)
    markerscoord = markers %*% t(v %*% solve(t(v) %*% v))
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
    geom_point(data = U, aes(x = .data$dim1, y = .data$dim2), colour = ocol) 
  # xlab(paste("Dimension", dims[1])) +
  # ylab(paste("Dimension", dims[2])) +
  # labs(title = "Type I plot")
  
  
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
  
  ######################################################
  # variable axes with ticks and markers for responses
  ######################################################
  plt = plt + geom_abline(intercept = 0, slope = V[,2]/V[,1], colour = ycol) +
    geom_label(data = MCy, aes(x = .data$dim1, y = .data$dim2, label = labs), 
               fill = ycol, fontface = "bold", color = "white", size = markersize)
  
  # geom_point(data = MCy, aes(x = .data$dim1, y = .data$dim2), shape = 15, colour = ycol) +
  # geom_text(data = MCy, aes(x = .data$dim1, y = .data$dim2, label = labs), nudge_y = -0.08, size = 1.5)
  
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

###################################################

plot.lpca2 <- function(object, dims = c(1,2), 
                       ycol = "darkgreen", xcol = "darkblue", ocol = "grey",
                       markersize = 2.5, labelsize = 3)
{
  # plots the results of a
  # logistic principal component analysis (X = NULL)
  # logistic reduced rank regression (X != NULL)
  # distance representation
  #
  # @param object an object of type lpca
  # @param dims which dimensions to visualize
  # @param ycol colour for representation of response variables
  # @param xcol colour for representation of predictor variables
  # @param ocol colour for representation of row objects
  #
  
  object = lpca2dist(object)
  
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
  rownamesVV<-rownames(VV)
  plt = ggplot() +
    geom_point(data = U, aes(x = .data$dim1, y = .data$dim2), colour = ocol) +
    geom_point(data = VV, aes(x = .data$dim1, y = .data$dim2), colour = ycol) +
    geom_text_repel(data = VV, aes(x = .data$dim1, y = .data$dim2, label = rownamesVV), family = "mono") 
  # xlab(paste("Dimension", dims[1])) +
  # ylab(paste("Dimension", dims[2])) +
  # labs(title = "Type D plot")
  
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
      geom_line(data = MCx1, aes(x = .data$dim1, y = .data$dim2, group = .data$varx), col = xcol) +
      geom_label(data = MCx2, aes(x = .data$dim1, y = .data$dim2, label = .data$labs), 
                 fill = xcol, fontface = "bold", color = "white", size = markersize) +
      geom_point(data = MCx3, aes(x = .data$dim1, y = .data$dim2), col = xcol, size = 4) +
      geom_text_repel(data = MCx3[-1, ], aes(x = .data$dim1, y = .data$dim2, label = labs,
                                             family = 'mono', fontface = 'bold'), size = 3)
  }
  
  ######################################################
  # variable labels for predictors
  ######################################################
  
  if(isx){
    BV = rbind(B[!dichotomous, ])
    PP = nrow(B[!dichotomous, ])
    CC = 0
    names = object$xnames[!dichotomous]
    
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
    
  }
  
  suppressWarnings(print(plt))
  return(plt)
}

plot.lpcah <- function(object, dims = c(1,2), pmarkers,
                       ycol = "darkgreen", xcol = "darkblue", ocol = "grey",
                       markersize = 2.5, labelsize = 3){
  # plots the results of a
  # logistic principal component analysis (X = NULL)
  # logistic reduced rank regression (X != NULL)
  #
  # @param object an object of type lpca
  # @param dims which dimensions to visualize
  # @param ycol colour for representation of response variables
  # @param xcol colour for representation of predictor variables
  # @param ocol colour for representation of row objects
  #
  
  ######################################################
  # retrieve information from object
  ######################################################
  object2 = lpca2dist(object)
  Y = object$Y
  
  U = as.data.frame(object$U[ , dims])
  N = nrow(U)
  colnames(U) = c("dim1", "dim2")
  
  V = object$V[ , dims]
  VV = as.data.frame(V)
  R = nrow(V)
  colnames(VV) = c("dim1", "dim2")
  rownames(VV) = object$ynames
  
  W = object2$V
  WW = cbind(W[seq(1,(nrow(W)-1),by = 2), ], W[seq(2, nrow(W), by = 2), ])
  WW = as.data.frame(WW)
  rownames(WW) = object$ynames
  colnames(WW) = c("dim1a", "dim2a", "dim1b", "dim2b")
  
  ######################################################
  # retrieve information for response variables variable axes
  ######################################################
  
  MCy <- data.frame(labs=character(),
                    vary = integer(),
                    dim1 = double(),
                    dim2 = double(), stringsAsFactors=FALSE)
  ll = 0
  
  if(!is.list(pmarkers)){
    pp = pmarkers
    pmarkers <- vector(mode = "list", length = R)
    for(r in 1:R){pmarkers[[r]] = pp}
  }
  
  for(r in 1:R){
    probs = pmarkers[[r]]
    lo = log(probs/(1 - probs))
    l.m = length(lo)
    markerlabs = paste(probs)
    
    markers = matrix(lo - object$m[r], length(lo), 1)
    v = matrix(V[r, dims], nrow = 2, ncol = 1)
    markerscoord = markers %*% t(v %*% solve(t(v) %*% v))
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
    geom_point(data = U, aes(x = .data$dim1, y = .data$dim2), colour = ocol, alpha = 0.5, size = 0.5) 
  # xlab(paste("Dimension", dims[1])) +
  # ylab(paste("Dimension", dims[2])) +
  # labs(title = "Hybrid triplot")
  
  
  ######################################################
  # variable axes with ticks and markers for predictors
  ######################################################
  if(isx){
    plt = plt + geom_abline(intercept = 0, slope = B[!dichotomous,2]/B[!dichotomous,1], colour = xcol, linetype = 3) +
      geom_line(data = MCx1, aes(x = .data$dim1, y = .data$dim2, group = .data$varx), col = xcol) +
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
    geom_abline(intercept = 0, slope = V[,2]/V[,1], colour = ycol, linetype = 3) +
    geom_label(data = MCy, aes(x = .data$dim1, y = .data$dim2, label = labs), 
               fill = ycol, fontface = "bold", color = "white", size = markersize)
  # geom_point(data = MCy, aes(x = .data$dim1, y = .data$dim2), size = 1, shape = 15, colour = ycol) +
  # geom_text(data = MCy, aes(x = .data$dim1, y = .data$dim2, label = labs), nudge_y = -0.08, size = 1.5)
  
  plt = plt + geom_segment(aes(x = .data$dim1a, y = .data$dim2a, xend = .data$dim1b, yend = .data$dim2b), colour = ycol, data = WW)
  
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

lpca2dist = function(object){
  # transformation of lpca object o distance represenations
  #
  K = -object$V/2
  B = object$B
  U = object$U
  a = object$m
  
  R = nrow(K)
  M = ncol(K)
  
  Ak = kronecker(diag(R), matrix(c(1,-1),2,1)) # for K - response variable discrimination
  Al = kronecker(diag(R), matrix(c(1,1),2,1)) # for L - response variable position
  
  L = matrix(0, R, M)
  
  # copy from melodic.R
  if(M == 1){L = a/K}
  if(M > 1){
    for(r in 1:R){
      L[r, ] = a[r] * K[r, ]/sum(K[r, ]^2)
    }
  }
  Vl = Al %*% L  # midpoints
  
  Vk = Ak %*% K # deviation from midpoints/discrimination
  Vl = Al %*% L  # midpoints/location
  V = Vl + Vk
  
  ynames = paste0(rep(object$ynames, each = 2), c(0,1))
  
  mldm = list(
    Y = object$Y,
    Xoriginal = object$Xoriginal,
    X = object$X,
    mx = object$mx,
    sdx = object$sdx,
    ynames = ynames,
    xnames = object$xnames,
    U = U,
    B = B,
    V = V,
    K = K,
    L = L,
    iter = object$iter,
    deviance = object$deviance
  )
  return(mldm)
}
