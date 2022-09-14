plot.lpca1 <- function(object, dims = c(1,2), ycol = "darkgreen", xcol = "lightskyblue", ocol = "grey")
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
  lo = log(seq(from = 0.1, to = 0.9, by = 0.1)/(1 - seq(from = 0.1, to = 0.9, by = 0.1)))
  l.m = length(lo)
  markerlabs = paste(seq(0.1, 0.9, by = 0.1))

  for(r in 1:R){
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
  } #isx

  ######################################################
  # plotting - objects
  ######################################################
  plt = ggplot() +
    geom_point(data = U, aes_string(x = 'dim1', y = 'dim2'), colour = ocol) +
    xlab(paste("Dimension", dims[1])) +
    ylab(paste("Dimension", dims[2]))

  ######################################################
  # variable axes with ticks and markers for predictors
  ######################################################
  if(isx){
    plt = plt + geom_abline(intercept = 0, slope = B[,2]/B[,1], colour = xcol, linetype = 3) +
      geom_line(data = MCx1, aes_string(x = 'dim1', y = 'dim2', group = 'varx'), col = xcol) +
      geom_point(data = MCx2, aes_string(x = 'dim1', y = 'dim2'), col = xcol) +
      geom_text(data = MCx2, aes_string(x = 'dim1', y = 'dim2', label = 'labs'), nudge_y = -0.08, size = 1.5)
  }

  ######################################################
  # variable axes with ticks and markers for responses
  ######################################################
  plt = plt + geom_abline(intercept = 0, slope = V[,2]/V[,1], colour = ycol) +
    geom_point(data = MCy, aes_string(x = 'dim1', y = 'dim2'), shape = 15, colour = ycol) +
    geom_text(data = MCy, aes_string(x = 'dim1', y = 'dim2', label = 'labs'), nudge_y = -0.08, size = 1.5)

  ######################################################
  # variable labels
  ######################################################
  a = ceiling(max(abs(c(ggplot_build(plt)$layout$panel_scales_x[[1]]$range$range, ggplot_build(plt)$layout$panel_scales_y[[1]]$range$range))))

  idx2 = apply(abs(V), 1, which.max)
  t2 = s2 = rep(NA,R)
  for(rr in 1:R){
    t2[rr] = (a * 1.1)/(abs(V[rr,idx2[rr]])) * V[rr,-idx2[rr]]
    s2[rr] = sign(V[rr,idx2[rr]])
  }
  CC2 = cbind(idx2, t2, s2)

  if(isx){
    idx1 = apply(abs(B), 1, which.max)
    t1 = s1 = rep(NA,P)
    for(pp in 1:P){
      t1[pp] = (a * 1.1)/(abs(B[pp,idx1[pp]])) * B[pp,-idx1[pp]]
      s1[pp] = sign(B[pp,idx1[pp]])
    }
    CC1 = cbind(idx1, t1, s1)
  }

  if(isx) CC = rbind(CC1, CC2) else CC = CC2
  if(isx) rownames(CC) = c(object$xnames, object$ynames) else rownames(CC) = object$ynames
  colnames(CC) = c("idx", "t", "s")

  bottom = which(CC[, "idx"] == 2 & CC[, "s"] == -1)
  top =  which(CC[, "idx"] == 2 & CC[, "s"] == 1)
  right = which(CC[, "idx"] == 1 & CC[, "s"] == 1)
  left = which(CC[, "idx"] == 1 & CC[, "s"] == -1)


  if(length(CC[top, "t"])==0){
    plt <- plt + scale_x_continuous(limits = c(-a,a), breaks = CC[bottom, "t"], labels = rownames(CC)[bottom],
                                    sec.axis = sec_axis(trans ~ ., breaks = 0, labels = ""))
  }else{
    plt <- plt + scale_x_continuous(limits = c(-a,a), breaks = CC[bottom, "t"], labels = rownames(CC)[bottom],
                                    sec.axis = sec_axis(trans ~ ., breaks = CC[top, "t"], labels = rownames(CC)[top]))
  }


  if(length(CC[right, "t"])==0){
    plt = plt + scale_y_continuous(limits = c(-a,a), breaks = CC[left, "t"], labels = rownames(CC)[left],
                                   sec.axis = sec_axis(trans ~ ., breaks = 0, labels = ""))
  }else{
    plt = plt + scale_y_continuous(limits = c(-a,a), breaks = CC[left, "t"], labels = rownames(CC)[left],
                                   sec.axis = sec_axis(trans ~ ., breaks = CC[right, "t"], labels = rownames(CC)[right]))
  }


  plt = plt + theme(axis.line = element_line(colour = "black"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank())
  plt = plt + coord_fixed()

  #suppressWarnings(print(plt))

  return(plt)
}





###################################################

plot.lpca2 <- function(object, dims = c(1,2), ycol = "darkgreen", xcol = "lightskyblue", ocol = "grey")
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
  # retrieve information for response variables variable axes
  ######################################################

  # MCy <- data.frame(labs=character(),
  #                   vary = integer(),
  #                   dim1 = double(),
  #                   dim2 = double(), stringsAsFactors=FALSE)
  # ll = 0
  # lo = log(seq(from = 0.1, to = 0.9, by = 0.1)/(1 - seq(from = 0.1, to = 0.9, by = 0.1)))
  # l.m = length(lo)
  # markerlabs = paste(seq(0.1, 0.9, by = 0.1))
  #
  # for(r in 1:R){
  #   markers = matrix(lo - object$m[r], length(lo), 1)
  #   v = matrix(V[r, dims], nrow = 2, ncol = 1)
  #   markerscoord = markers %*% t(v %*% solve(t(v) %*% v))
  #   MCy[(ll + 1): (ll + l.m), 1] = markerlabs
  #   MCy[(ll + 1): (ll + l.m), 2] = r
  #   MCy[(ll + 1): (ll + l.m), 3:4] = markerscoord
  #   ll = ll + l.m
  # }
  #
  ######################################################
  # retrieve information for predictor variables variable axes
  ######################################################
  X = object$X
  isx = !is.null(X)
  if(isx){
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
  } #isx

  ######################################################
  # plotting - objects
  ######################################################
  plt = ggplot() +
    geom_point(data = U, aes_string(x = 'dim1', y = 'dim2'), colour = ocol) +
    # geom_point(data = U, aes(x = dim1, y = dim2), colour = ocol, shape = ".") +
    geom_point(data = VV, aes_string(x = 'dim1', y = 'dim2'), colour = ycol) +
    geom_text_repel(data = VV, aes_string(x = 'dim1', y = 'dim2', label = 'rownames(VV)'), family = "mono") +
    #geom_text(data = VV, aes(x = dim1, y = dim2, label = rownames(VV)), family = "mono") +
    xlab(paste("Dimension", dims[1])) +
    ylab(paste("Dimension", dims[2]))

  ######################################################
  # variable axes with ticks and markers for predictors
  ######################################################
  if(isx){
    plt = plt + geom_abline(intercept = 0, slope = B[,2]/B[,1], colour = xcol, linetype = 3) +
      geom_line(data = MCx1, aes_string(x = 'dim1', y = 'dim2', group = 'varx'), col = xcol) +
      geom_point(data = MCx2, aes_string(x = 'dim1', y = 'dim2'), col = xcol) +
      geom_text(data = MCx2, aes_string(x = 'dim1', y = 'dim2', label = 'labs'), nudge_y = -0.08, size = 1.5)
  }

  ######################################################
  # variable axes with ticks and markers for responses
  ######################################################
  # plt = plt + geom_abline(intercept = 0, slope = V[,2]/V[,1], colour = ycol) +
  #   geom_point(data = MCy, aes(x = dim1, y = dim2), shape = 15, colour = ycol) +
  #   geom_text(data = MCy, aes(x = dim1, y = dim2, label = labs), nudge_y = -0.08, size = 1.5)

  ######################################################
  # variable labels
  ######################################################
  a = ceiling(max(abs(c(ggplot_build(plt)$layout$panel_scales_x[[1]]$range$range, ggplot_build(plt)$layout$panel_scales_y[[1]]$range$range))))

  # idx2 = apply(abs(V), 1, which.max)
  # t2 = s2 = rep(NA,R)
  # for(rr in 1:R){
  #   t2[rr] = (a * 1.1)/(abs(V[rr,idx2[rr]])) * V[rr,-idx2[rr]]
  #   s2[rr] = sign(V[rr,idx2[rr]])
  # }
  # CC2 = cbind(idx2, t2, s2)

  if(isx){
    idx1 = apply(abs(B), 1, which.max)
    t1 = s1 = rep(NA,P)
    for(pp in 1:P){
      t1[pp] = (a * 1.1)/(abs(B[pp,idx1[pp]])) * B[pp,-idx1[pp]]
      s1[pp] = sign(B[pp,idx1[pp]])
    }
    CC = cbind(idx1, t1, s1)
    rownames(CC) = object$xnames
    colnames(CC) = c("idx", "t", "s")

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


      }

  plt = plt + theme(axis.line = element_line(colour = "black"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank())
  plt = plt + coord_fixed()

  #suppressWarnings(print(plt))

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
