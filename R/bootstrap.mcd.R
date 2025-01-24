#'  Bootstrap procedure for Multonimal Canonical Decomposition Model
#'
#' @param object An output object from mcd1 or mcd2
#' @param Bsamples Number of Bootstrap samples to take
#'
#' @return BBdf Bootstrap estimates of B
#' @return BVdf Bootstrap estimates of V
#'
#' @examples
#' \dontrun{
#' data(dataExample_lpca)
#' Y = as.matrix(dataExample_lpca[, 1:8])
#' X = as.matrix(dataExample_lpca[, 9:13])
#' # supervised
#' output = mcd1(X, Y, S = 2, ord.z = 2)
#' #' boot.output = bootstrap.mcd(output, Bsamples = 100)
#' plot(boot.output)
#' }
#'
#' @export

bootstrap.mcd = function(object, Bsamples = 1000){
  if(object$type == "mcd1"){
    bres = bootstrap.mcd1(object, Bsamples)
  }
  else if (object$type == "mcd2"){
    bres = bootstrap.mcd2(object, Bsamples)
  }
  return(bres)
}

# for multivariate binary

bootstrap.mcd1 = function(object, Bsamples = 1000){
  # performs a bootstrap for a model fitted with the mcd() function

  X = object$X
  Y = object$Y
  Z = object$Z
  W = object$W
  qrz = object$qrz
  G = object$G
  A = object$A

  N = nrow(X)
  S = ncol(object$Bx)
  R = nrow(object$Bz)
  P = ncol(X)
  TT = ncol(Z)
  Q = ncol(W)

  # balanced bootstrap scheme
  f = matrix(1:N, N, Bsamples)
  ff = matrix(f,prod(dim(f)),1)
  fff = sample(ff)
  f = matrix(fff, N, Bsamples)

  # starting values for bootstrap analyses
  start = list(m = object$m,
               bm = object$bm,
               Bx = object$Bx,
               Bz = object$Bz)


  # create empty matrices for bootstrap estimates
  BBx = matrix(NA, P*S, Bsamples)
  BBz = matrix(NA, TT*S, Bsamples)
  BA = matrix(NA, P * TT, Bsamples)
  Bm = matrix(NA, length(object$m), Bsamples)
  Bbm = matrix(NA, length(object$bm), Bsamples)
  BBxdf = matrix(NA, P*Bsamples, (S + 2))
  BBzdf = matrix(NA, TT*Bsamples, (S + 2))

  for(b in 1:Bsamples){
    cat("This is analysis", b, "from a total of", Bsamples, "Bootstraps", "\n")
    obs <- f[ , b]
    bres = bmcd1(X[obs, ], Y[obs, ], Z, W, G[obs, ], A, qrz, start)

    #
    Bm[, b] = bres$m
    Bbm[, b] = bres$bm

    BBx[, b] = matrix(bres$Bx, ncol = 1)
    BBz[, b] = matrix(bres$Bz, ncol = 1)
    BA[ , b] = matrix((bres$Bx %*% t(bres$Bz)), ncol = 1)

    BBxdf[((b-1)*P + 1):(b*P), 1] = b
    BBxdf[((b-1)*P + 1):(b*P), 2] = 1:P
    BBxdf[((b-1)*P + 1):(b*P), 3:(S+2)] = bres$Bx

    BBzdf[((b-1)*TT + 1):(b*TT), 1] = b
    BBzdf[((b-1)*TT + 1):(b*TT), 2] = 1:TT
    BBzdf[((b-1)*TT + 1):(b*TT), 3:(S+2)] = bres$Bz

  }

  se.A =  matrix(apply(X = BA, MARGIN = 1, FUN = "sd"), ncol = TT)
  rownames(se.A) = object$xnames
  colnames(se.A) = object$znames

  BBxdf = as.data.frame(BBxdf)
  colnames(BBxdf) = c("Bootstrap", "Predictor", paste0("dim", 1:S))
  BBxdf$Predictor = factor(BBxdf$Predictor, levels = 1:P, labels = object$xnames)

  BBzdf = as.data.frame(BBzdf)
  colnames(BBzdf) = c("Bootstrap", "Response", paste0("dim", 1:S))
  BBzdf$Response = factor(BBzdf$Response, levels = 1:TT, labels = object$znames)

  b.output = list(
    mcdobj = object,
    BBx = BBx,
    BBz = BBz,
    BA = BA,
    se.A = se.A, # standard deviations of bootstrap coefficients
    BBxdf = BBxdf,
    BBzdf = BBzdf,
    Bbm = Bbm,
    Bm = Bm
  )
  class(b.output) = "bootstrap"
  return(b.output)
}

bmcd1 = function(X, Y, Z, W, G, A, qrz, start, trace = FALSE, maxiter = 65536, dcrit = 1e-6){
  # bootstrap version for multinomial canonical decomposition

  n = nrow(X)
  ones.n = matrix(1,n,1)
  P = ncol(X)
  R = ncol(Y)
  S = ncol(start$Bx)
  Q = ncol(Z)

  if(P == 1){
    iRx = 1/sqrt(t(X) %*% X)
  }
  else{
    eig.out = eigen(t(X) %*% X)
    iRx = eig.out$vectors %*% diag(1/sqrt(eig.out$values)) %*% t(eig.out$vectors)
  }

  if(ncol(Z) == 1){
    iRz = 1/sqrt(t(Z) %*% Z)
  }
  else{
    eig.out = eigen(t(Z) %*% Z)
    iRz = eig.out$vectors %*% diag(1/sqrt(eig.out$values)) %*% t(eig.out$vectors)
  }

  # initialization
  m = start$m
  bm = start$bm
  #
  Bx = start$Bx
  Bz = start$Bz
  theta = ones.n %*% t(m) + X %*% Bx %*% t(Z %*% Bz)
  Ghat = exp(theta) / rowSums(exp(theta))
  dev.old = -2 * sum(log(Ghat[G == 1]))

  # iteration
  iter = 0; dif = 1
  while(dif > dcrit){
    iter = iter + 1
    # update m
    ZZ = theta + 4 * (G - Ghat)
    m = colMeans(ZZ - X %*% Bx %*% t(Z %*% Bz))
    # bm = iWWW %*% m
    # m = W %*% bm
    bm = solve.qr(qrz, m)
    m = qr.fitted(qrz, m)

    # update Bx (U) and Bz (V)
    theta = ones.n %*% t(m) + X %*% Bx %*% t(Z %*% Bz)
    Ghat = exp(theta) / rowSums(exp(theta))
    ZZ = theta + 4 * (G - Ghat) - ones.n %*% t(m)
    udv = svd(iRx %*% t(X) %*% ZZ %*% Z %*% iRz)
    Bx = iRx %*% matrix(udv$u[, 1:S], P, S) %*% diag(udv$d[1:S], nrow = S, ncol = S)
    Bz = iRz %*% matrix(udv$v[, 1:S], Q, S)

    # deviance
    theta = ones.n %*% t(m) + X %*% Bx %*% t(Z %*% Bz)
    Ghat = exp(theta) / rowSums(exp(theta))
    dev.new = -2 * sum(log(Ghat[G == 1]))

    # convergence
    dif = 2 * (dev.old - dev.new)/ ((dev.old + dev.new))
    if ( trace ) cat(iter, dev.old, dev.new, dif, "\n")
    if ( dif < dcrit ) break
    if ( iter == maxiter ) warning("Maximum number of iterations reached - not converged (yet)")
    dev.old = dev.new
  } # end iteration

  maxs = apply(start$Bx, 2, which.max)
  for(s in 1:S){
    if(Bx[maxs[s], s] < 0){
      Bx[, s] = -Bx[, s]
      Bz[, s] = -Bz[, s]
    }
  }

  BB = Bx %*% t(Bz)

  # name giving
  rownames(Bz) = colnames(Z)
  rownames(Bx) = colnames(X)
  rownames(BB) = colnames(X)
  colnames(BB) = colnames(Z)
  rownames(bm) = colnames(W)

  npar = length(bm) + (P + Q - S) * S
  # create output object
  results = list(
    m = m,
    bm = bm,
    Bx = Bx,
    Bz = Bz
  )
  return(results)
}

# for multinomial

bootstrap.mcd2 = function(object, Bsamples = 1000){
  # performs a bootstrap for a model fitted with the mcd() function

  X = object$X
  Z = object$Z
  G = object$G

  N = nrow(X)
  S = ncol(object$Bx)
  R = nrow(object$Bz)
  P = ncol(X)
  TT = ncol(Z)

  # balanced bootstrap scheme
  f = matrix(1:N, N, Bsamples)
  ff = matrix(f,prod(dim(f)),1)
  fff = sample(ff)
  f = matrix(fff, N, Bsamples)

  # starting values for bootstrap analyses
  start = list(m = object$m,
               bm = object$bm,
               Bx = object$Bx,
               Bz = object$Bz)


  # create empty matrices for bootstrap estimates
  BBx = matrix(NA, P*S, Bsamples)
  BBz = matrix(NA, TT*S, Bsamples)
  BA = matrix(NA, P * TT, Bsamples)
  Bm = matrix(NA, length(object$m), Bsamples)
  # Bbm = matrix(NA, length(object$bm), Bsamples)
  BBxdf = matrix(NA, P*Bsamples, (S + 2))
  BBzdf = matrix(NA, TT*Bsamples, (S + 2))

  for(b in 1:Bsamples){
    cat("This is analysis", b, "from a total of", Bsamples, "Bootstraps", "\n")
    obs <- f[ , b]
    bres = bmcd2(X[obs, ], G[obs, ], Z, start)

    #
    Bm[, b] = bres$m
    # Bbm[, b] = bres$bm

    BBx[, b] = matrix(bres$Bx, ncol = 1)
    BBz[, b] = matrix(bres$Bz, ncol = 1)
    BA[ , b] = matrix((bres$Bx %*% t(bres$Bz)), ncol = 1)

    BBxdf[((b-1)*P + 1):(b*P), 1] = b
    BBxdf[((b-1)*P + 1):(b*P), 2] = 1:P
    BBxdf[((b-1)*P + 1):(b*P), 3:(S+2)] = bres$Bx

    BBzdf[((b-1)*TT + 1):(b*TT), 1] = b
    BBzdf[((b-1)*TT + 1):(b*TT), 2] = 1:TT
    BBzdf[((b-1)*TT + 1):(b*TT), 3:(S+2)] = bres$Bz

  }

  se.A =  matrix(apply(BA, 1, "sd"), ncol = TT)
  rownames(se.A) = object$xnames
  colnames(se.A) = object$znames

  BBxdf = as.data.frame(BBxdf)
  colnames(BBxdf) = c("Bootstrap", "Predictor", paste0("dim", 1:S))
  BBxdf$Predictor = factor(BBxdf$Predictor, levels = 1:P, labels = object$xnames)

  BBzdf = as.data.frame(BBzdf)
  colnames(BBzdf) = c("Bootstrap", "Response", paste0("dim", 1:S))
  BBzdf$Response = factor(BBzdf$Response, levels = 1:TT, labels = object$znames)

  b.output = list(
    mcdobj = object,
    BBx = BBx,
    BBz = BBz,
    BA = BA,
    se.A = se.A, # standard deviations of bootstrap coefficients
    BBxdf = BBxdf,
    BBzdf = BBzdf,
    # Bbm = Bbm,
    Bm = Bm
  )
  class(b.output) = "bootstrap"
  return(b.output)
}

bmcd2 = function(X, G, Z, start, trace = FALSE, maxiter = 65536, dcrit = 1e-6){
  # bootstrap version for multinomial canonical decomposition 2

  n = nrow(X)
  ones.n = matrix(1,n,1)
  P = ncol(X)
  S = ncol(start$Bx)
  Q = ncol(Z)

  if(P == 1){
    iRx = 1/sqrt(t(X) %*% X)
  }
  else{
    eig.out = eigen(t(X) %*% X)
    iRx = eig.out$vectors %*% diag(1/sqrt(eig.out$values)) %*% t(eig.out$vectors)
  }

  if(ncol(Z) == 1){
    iRz = 1/sqrt(t(Z) %*% Z)
  }
  else{
    eig.out = eigen(t(Z) %*% Z)
    iRz = eig.out$vectors %*% diag(1/sqrt(eig.out$values)) %*% t(eig.out$vectors)
  }

  # initialization
  m = start$m
  #
  Bx = start$Bx
  Bz = start$Bz
  theta = ones.n %*% t(m) + X %*% Bx %*% t(Z %*% Bz)
  Ghat = exp(theta) / rowSums(exp(theta))
  dev.old = -2 * sum(log(Ghat[G == 1]))

  # iteration
  iter = 0; dif = 1
  while(dif > dcrit){
    iter = iter + 1
    # update m
    ZZ = theta + 4 * (G - Ghat)
    m = colMeans(ZZ - X %*% Bx %*% t(Z %*% Bz))

    # update Bx (U) and Bz (V)
    theta = ones.n %*% t(m) + X %*% Bx %*% t(Z %*% Bz)
    Ghat = exp(theta) / rowSums(exp(theta))
    ZZ = theta + 4 * (G - Ghat) - ones.n %*% t(m)
    udv = svd(iRx %*% t(X) %*% ZZ %*% Z %*% iRz)
    Bx = iRx %*% matrix(udv$u[, 1:S], P, S) %*% diag(udv$d[1:S], nrow = S, ncol = S)
    Bz = iRz %*% matrix(udv$v[, 1:S], Q, S)

    # deviance
    theta = ones.n %*% t(m) + X %*% Bx %*% t(Z %*% Bz)
    Ghat = exp(theta) / rowSums(exp(theta))
    dev.new = -2 * sum(log(Ghat[G == 1]))

    # convergence
    dif = 2 * (dev.old - dev.new)/ ((dev.old + dev.new))
    if ( trace ) cat(iter, dev.old, dev.new, dif, "\n")
    if ( dif < dcrit ) break
    if ( iter == maxiter ) warning("Maximum number of iterations reached - not converged (yet)")
    dev.old = dev.new
  } # end iteration

  maxs = apply(start$Bx, 2, which.max)
  for(s in 1:S){
    if(Bx[maxs[s], s] < 0){
      Bx[, s] = -Bx[, s]
      Bz[, s] = -Bz[, s]
    }
  }

  BB = Bx %*% t(Bz)

  # name giving
  rownames(Bz) = colnames(Z)
  rownames(Bx) = colnames(X)
  rownames(BB) = colnames(X)
  colnames(BB) = colnames(Z)

  npar = length(m) + (P + Q - S) * S
  # create output object
  results = list(
    m = m,
    Bx = Bx,
    Bz = Bz
  )
  return(results)
}

