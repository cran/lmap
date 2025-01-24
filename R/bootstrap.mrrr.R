#'  Bootstrap procedure for Multinomial Reduced Rank Model
#'
#' @param object An output object from mrrr
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
#' output = mrrr(y = y, X = X, S = 2)
#' boot.output = bootstrap.mrrr(output, Bsamples = 100)
#' plot(boot.output)
#' }
#'
#' @export

bootstrap.mrrr = function(object, Bsamples = 1000, myseed = 1){

  N = nrow(object$G)
  X = object$X
  G = object$G
  y = object$y
  bstart = list(m = object$m, B = object$B, V = object$V)
  P = ncol(X)
  R = nrow(object$V)
  S = ncol(object$V)
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
    # determine bootstrap sample and refit model
    obs <- f[ , b]
    b.out = mrrr(y = y[obs], X = X[obs, ], S = S, trace = F, start = bstart)

    # out of sample prediction - quality measure is average deviance
    Ghat = predict.mrrr(b.out, X[-obs, ])
    sdev.oos[b] = -2 * mean(log(Ghat[G[-obs, ] == 1]))

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
