#'  Bootstrap procedure for Multinomial Restricted Unfolding
#'
#' @param object An output object from mru
#' @param Bsamples Number of Bootstrap samples to take
#' @param myseed A seed number to make the bootstrap reproducible
#'
#' @return BBdf Bootstrap estimates of B
#' @return BVdf Bootstrap estimates of V
#'
#' @examples
#' \dontrun{
#' data(dataExample_mru)
#' y = as.matrix(dataExample_mru[ , 1])
#' X = as.matrix(dataExample_mru[ , 2:6])
#' output = mru(y = y, X = X, S = 2)
#' boot.output = bootstrap.mru(output, Bsamples = 100)
#' plot(boot.output)
#' }
#'
#' @export

bootstrap.mru = function(object, Bsamples = 1000, myseed = 1){

  G = object$G
  N = nrow(object$G)
  X = object$X
  VV = object$V
  bstart = list(B = object$Bx, V = object$V)
  P = ncol(X)
  R = nrow(object$V)
  S = ncol(object$V)
  y = object$y
  object$B = object$Bx # trick for plot.bootstrap function
  
  set.seed(myseed)

  # balanced bootstrap scheme
  if(is.matrix(Bsamples)){f = Bsamples; Bsamples = ncol(f)}
  # else create balanced bootstrap samples
  else{
    f = matrix(1:N, N, Bsamples)
    ff = matrix(f,prod(dim(f)),1)
    fff = sample(ff)
    f = matrix(fff, N, Bsamples)
  }


  # create empty matrices for bootstrap estimates
  BB = matrix(NA, P*S, Bsamples)
  BV = matrix(NA, R*S, Bsamples)
  BBdf = matrix(NA, P*Bsamples, (S + 2))
  BVdf = matrix(NA, R*Bsamples, (S + 2))
  sdev.oos = rep(NA, Bsamples)

  # bootstrap analysis
  for(b in 1:Bsamples){
    cat("This is analysis", b, "from a total of", Bsamples, "Bootstraps", "\n")
    obs <- f[ , b]
    b.out = mru(y = G[obs, ], X = X[obs, ], S = S, start = bstart)

    # out of sample prediction - quality measure is average deviance
    Ghat = predict.mru(b.out, X[-obs, ])$Yhat
    sdev.oos[b] = -2 * mean(log(Ghat[G[-obs, ] == 1]))

    # procrustes rotatie op V
    pq = svd(t(VV) %*% b.out$V)
    TT = pq$v %*% t(pq$u)

    # rotatie van schattingen
    b.out$V = b.out$V %*% TT
    b.out$Bx = b.out$Bx %*% TT
    b.out$U = b.out$U %*% TT

    # collect boostrap estimates
    BV[, b] = matrix(b.out$V, ncol = 1)
    BVdf[((b-1)*R + 1):(b*R), 1] = b
    BVdf[((b-1)*R + 1):(b*R), 2] = 1:R
    BVdf[((b-1)*R + 1):(b*R), 3:(S+2)] = b.out$V

    BB[, b] = matrix(b.out$Bx, ncol = 1)
    BBdf[((b-1)*P + 1):(b*P), 1] = b
    BBdf[((b-1)*P + 1):(b*P), 2] = 1:P
    BBdf[((b-1)*P + 1):(b*P), 3:(S+2)] = b.out$Bx
  }

  BBdf = as.data.frame(BBdf)
  colnames(BBdf) = c("Bootstrap", "Predictor", paste0("dim", 1:S))
  BBdf$Predictor = factor(BBdf$Predictor, levels = 1:P, labels = object$xnames)

  BVdf = as.data.frame(BVdf)
  colnames(BVdf) = c("Bootstrap", "Response", paste0("dim", 1:S))
  BVdf$Response = factor(BVdf$Response, levels = 1:R, labels = object$ynames)

  output = list(
    obj = object,
    BB = BB,
    BBdf = BBdf,
    BV = BV,
    BVdf = BVdf,
    sdev.oos = sdev.oos
  )

  class(output) = "bootstrap"
  return(output)
}
