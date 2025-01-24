#'  Bootstrap procedure for Logistic (Restricted) MDU
#'
#' @param object An output object from lmdu
#' @param Bsamples Number of Bootstrap samples to take
#' @param myseed A seed number to make the bootstrap reproducible
#'
#' @return BBdf Bootstrap estimates of B
#' @return BVdf Bootstrap estimates of V
#'
#' @examples
#' \dontrun{
#' data(dataExample_lmdu)
#' Y = as.matrix(dataExample_lmdu[ , 1:8])
#' X = as.matrix(dataExample_lmdu[ , 9:13])
#' output2 = lmdu(Y = Y, X = X, S = 2)
#' boot.output = bootstrap.lmdu(output2, Bsamples = 100)
#' plot(boot.output)
#' }
#'
#' @export

bootstrap.lmdu = function(object, Bsamples = 1000, myseed = 1){

  N = nrow(object$Y)
  Y = object$Y
  VV = object$V
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

  # start bootstrap

  if(is.null(object$X)){
    bstart = list(U = object$U, V = object$V, m = object$m)
    R = nrow(object$V)
    S = ncol(object$V)

    # create empty matrices for bootstrap estimates
    BV = matrix(NA, R * S, Bsamples)
    BVdf = matrix(NA, R * Bsamples, (S + 2))

    BB = BBdf = sdev.oos = NULL

    # bootstrap analysis
    for(b in 1:Bsamples){
      cat("This is analysis", b, "from a total of", Bsamples, "Bootstraps", "\n")
      obs <- f[ , b]
      b.out = lmdu(Y = Y[obs, ], S = S, start = bstart)
      # procrustes rotatie op V
      pq = svd(t(VV) %*% b.out$V)
      TT = pq$v %*% t(pq$u)

      # rotatie van schattingen
      b.out$V = b.out$V %*% TT
      b.out$U = b.out$U %*% TT

      # collect boostrap estimates
      BV[, b] = matrix(b.out$V, ncol = 1)
      BVdf[((b-1)*R + 1):(b*R), 1] = b
      BVdf[((b-1)*R + 1):(b*R), 2] = 1:R
      BVdf[((b-1)*R + 1):(b*R), 3:(S+2)] = b.out$V
    }
  }

  if(!is.null(object$X)){
    X = object$X
    bstart = list(B = object$B, V = object$V, m = object$m)
    P = ncol(X)
    R = nrow(object$V)
    S = ncol(object$V)
    Q = 2*Y - 1

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
      b.out = lmdu(Y = Y[obs, ], X = X[obs, ], S = S, start = bstart)

      # out of sample prediction - quality measure is average deviance
      theta = predict.lmdu(b.out, X[-obs, ])$theta
      sdev.oos[b] <- -2 * mean(log(plogis(Q[-obs, ] * theta)))

      # procrustes rotatie op V
      pq = svd(t(VV) %*% b.out$V)
      TT = pq$v %*% t(pq$u)

      # rotatie van schattingen
      b.out$V = b.out$V %*% TT
      b.out$B = b.out$B %*% TT
      b.out$U = b.out$U %*% TT


      # collect boostrap estimates
      BV[, b] = matrix(b.out$V, ncol = 1)
      BVdf[((b-1)*R + 1):(b*R), 1] = b
      BVdf[((b-1)*R + 1):(b*R), 2] = 1:R
      BVdf[((b-1)*R + 1):(b*R), 3:(S+2)] = b.out$V

      BB[, b] = matrix(b.out$B, ncol = 1)
      BBdf[((b-1)*P + 1):(b*P), 1] = b
      BBdf[((b-1)*P + 1):(b*P), 2] = 1:P
      BBdf[((b-1)*P + 1):(b*P), 3:(S+2)] = b.out$B
    }
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
