#' Multinomial Logistic Regression
#'
#' The function mlr performs multinomial logistic regression
#' for a nominal response variable and a set of predictor variables.
#' It uses an MM algorithm
#'
#' @param y An N vector of the responses (categorical).
#' @param X An N by P matrix with predictor variables
#' @param base The category that should be used as baseline. Can be NULL, in which case the colmeans are equal to zero. Can also be "largest", in which case the
#' @param maxiter maximum number of iterations
#' @param dcrit convergence criterion
#' @return Xoriginal Matrix X from input
#' @return X Scaled X matrix
#' @return G class indicator matrix
#' @return ynames class names of response variable
#' @return xnames variable names of the predictors
#' @return mx means of the predictor variables
#' @return sdx standard deviations of the predictor variables
#' @return A matrix with regression coefficients
#' @return iter number of iterations
#' @return deviance value of the deviance at convergence
#'
#' @examples
#' \dontrun{
#' data(dataExample_mru)
#' y = as.matrix(dataExample_mru[ , 1])
#' X = as.matrix(dataExample_mru[ , 2:6])
#' output = mlr(y = y, X = X, base = 1)
#' }
#'
#' @importFrom nnet class.ind
#'
#' @export

mlr = function(y, X, base = "largest", maxiter = 65536, dcrit = 1e-6){

  # preparation
  cal = match.call()

  trace = FALSE
  n = nrow(X)
  P = ncol(X)
  G = class.ind(y)
  C = ncol(G)

  if(!is.null(base)){if(base == "largest"){base = which.max(colSums(G))}}

  ones.C = matrix(1,C,1)
  Jc = diag(C) - 1/C
  ones.n = matrix(1,n,1)
  Jn = diag(n) - 1/n

  # scaling of predictor matrices
  Xoriginal = X
  outx = procx(X)
  X = outx$X
  mx = outx$mx
  sdx = outx$sdx
  P = ncol(X)

  xnames = colnames(X)
  ynames = colnames(G)

  # eig.out = eigen(t(X) %*% X)
  iXXX = solve(t(X) %*% X) %*% t(X)
  # iRx = eig.out$vectors %*% diag(1/sqrt(eig.out$values)) %*% t(eig.out$vectors)
  # iRxX = iRx %*% t(X)
  #
  # iGGG = solve(t(G) %*% G) %*% t(G)
  # A = t(t(X) %*% G %*% solve(t(G) %*% G))
  # iAAA = A %*% solve(t(A) %*% A)
  # eig.out = eigen(t(A) %*% A)
  # iRa = eig.out$vectors %*% diag(1/sqrt(eig.out$values)) %*% t(eig.out$vectors)
  # iRaA = iRa %*% t(A)

  # initialization
  m = colMeans(G)
  theta = ones.n %*% t(m)
  Ghat = exp(theta) / rowSums(exp(theta))
  Z = theta + 4 * (G - Ghat)
  A = iXXX %*% Z

  # # volgende kan wel weg binnenkort
  # if(random){
  #   B = matrix(rnorm(prod(dim(B))), nrow(B), ncol(B))
  #   V = matrix(rnorm(prod(dim(V))), nrow(V), ncol(V))
  # }

  theta = ones.n %*% t(m) + X %*% A
  Ghat = exp(theta) / rowSums(exp(theta))
  dev.old = -2 * sum(log(Ghat[G == 1]))

  # iteration
  iter = 0; dif = 1
  while(iter < maxiter){
    iter = iter + 1
    # update m
    Z = theta + 4 * (G - Ghat)

    # update m
    m = colMeans((Z - X %*% A) %*% Jc)

    # update B and V
    theta = ones.n %*% t(m) + X %*% A
    Ghat = exp(theta) / rowSums(exp(theta))
    Z = (theta + 4 * (G - Ghat)) %*% Jc
    A = iXXX %*% Z

    # udv = svd(iRxX %*% scale(Z, center = m, scale = FALSE) %*% Jc)
    # B = iRx %*% matrix(udv$u[, 1:S], P, S) %*% diag(udv$d[1:S], nrow = S, ncol = S)
    # V = matrix(udv$v[, 1:S], C, S)

    # deviance
    theta = ones.n %*% t(m) + X %*% A
    Ghat = exp(theta) / rowSums(exp(theta))
    dev.new = -2 * sum(log(Ghat[G == 1]))

    # convergence
    dif = 2 * (dev.old - dev.new)/(dev.old + dev.new)
    if ( trace ) cat(iter, dev.old, dev.new, dif, "\n")
    if ( dif < dcrit ) break
    if ( iter == maxiter ) warning("Maximum number of iterations reached - not converged (yet)")
    dev.old = dev.new
  } # end iteration

  # doe iets met baseline
  if(!is.null(base)){
    m = m - m[base]
    A = A - outer(A[ , base], rep(1, C))
  }
  rownames(A) = xnames
  colnames(A) = ynames

  theta = ones.n %*% t(m) + X %*% A
  Ghat = exp(theta) / rowSums(exp(theta))
  dev.new = -2 * sum(log(Ghat[G == 1]))

  npar = C - 1 + P * (C - 1)

  # create output object
  results = list(
    y = y,
    Xoriginal = Xoriginal,
    X = X,
    mx = mx,
    sdx = sdx,
    xnames = xnames,
    ynames = ynames,
    G = G,
    m = m,
    A = A,
    Ghat = Ghat,
    deviance = dev.new,
    npar = npar,
    AIC = dev.new + 2 * npar,
    BIC = dev.new + log(n) * npar,
    iter = iter
  )
  class(results) = "mlr"
  return(results)
}

