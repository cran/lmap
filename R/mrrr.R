#' Multinomial Reduced Rank Regression
#'
#' The function mrrr performs multinomial reduced rank regression for a nominal response
#' variable and a set of predictor variables.
#'
#' @param y An N vector of the responses (categorical).
#' @param X An N by P matrix with predictor variables
#' @param S Positive number indicating the dimensionality of teh solution
#' @param trace Boolean indicating whether a trace of the algorithm should be printed on the console.
#' @param maxiter maximum number of iterations
#' @param dcrit convergence criterion
#' @param start start values. If start=NULL, the algorithm computes the start values.
#' @return Xoriginal Matrix X from input
#' @return X Scaled X matrix
#' @return G class indicator matrix
#' @return ynames class names of response classes
#' @return xnames variable names of the predictors
#' @return mx means of the predictor variables
#' @return sdx standard deviations of the predictor variables
#' @return U coordinate matrix of row objects
#' @return B matrix with regression coefficients
#' @return V Class coordinate matrix
#' @return iters number of iterations
#' @return deviance value of the deviance at convergence
#'
#' @examples
#' \dontrun{
#' data(dataExample_mru)
#' y = as.matrix(dataExample_mru[ , 1])
#' X = as.matrix(dataExample_mru[ , 2:6])
#' output = mrrr(y = y, X = X, S = 2)
#' }
#'
#' @importFrom nnet class.ind
#'
#' @export

mrrr = function(y, X, S = 2, trace = FALSE, maxiter = 65536, dcrit = 1e-6, start = NULL){

  # preparation
  cal = match.call()

  n = nrow(X)
  P = ncol(X)
  G = class.ind(y)
  C = ncol(G)

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
  
  eig.out = eigen(t(X) %*% X)
  iXXX = solve(t(X) %*% X) %*% t(X)
  iRx = eig.out$vectors %*% diag(1/sqrt(eig.out$values)) %*% t(eig.out$vectors)
  iRxX = iRx %*% t(X)


  # initialization
  if(is.null(start)){
    m = colMeans(G)
    udv = svd(1/sqrt(n) * iRx %*% t(X) %*% G)
    B = iRx %*% matrix(udv$u[, 1:S], P, S) * sqrt(n)
    V = matrix(udv$v[, 1:S], C, S) %*% diag(udv$d[1:S], nrow = S, ncol = S) / sqrt(n)
  }
  else{
    m = start$m
    B = start$B
    V = start$V
  }

  theta = ones.n %*% t(m) + X %*% B %*% t(V)
  Ghat = exp(theta) / rowSums(exp(theta))
  dev.old = -2 * sum(log(Ghat[G == 1]))

  # iteration
  iter = 0; dif = 1
  while(iter < maxiter){
    iter = iter + 1
    # update m
    Z = theta + 4 * (G - Ghat)
    m = colMeans((Z - X %*% B %*% t(V)) %*% Jc)

    # update B and V
    theta = ones.n %*% t(m) + X %*% B %*% t(V)
    Ghat = exp(theta) / rowSums(exp(theta))
    Z = theta + 4 * (G - Ghat)
    udv = svd(iRxX %*% scale(Z, center = m, scale = FALSE) %*% Jc)
    B = iRx %*% matrix(udv$u[, 1:S], P, S) %*% diag(udv$d[1:S], nrow = S, ncol = S)
    V = matrix(udv$v[, 1:S], C, S)

    # deviance
    theta = ones.n %*% t(m) + X %*% B %*% t(V)
    Ghat = exp(theta) / rowSums(exp(theta))
    dev.new = -2 * sum(log(Ghat[G == 1]))

    # convergence
    dif = 2 * (dev.old - dev.new)/(dev.old + dev.new)
    if ( trace ) cat(iter, dev.old, dev.new, dif, "\n")
    if ( dif < dcrit ) break
    if ( iter == maxiter ) warning("Maximum number of iterations reached - not converged (yet)")
    dev.old = dev.new
  } # end iteration

  npar = C-1 + (P + C - S)*S
  rownames(B) = colnames(X)
  rownames(V) = colnames(G)

  # create output object
  results = list(
    call = cal,
    Xoriginal = Xoriginal,
    X = X,
    y = y,
    mx = mx,
    sdx = sdx,
    xnames = colnames(X),
    ynames = colnames(G),
    G = G,
    m = m,
    B = B,
    U = X %*% B,
    V = V,
    svd = udv,
    Ghat = Ghat,
    deviance = dev.new,
    npar = npar,
    AIC = dev.new + 2 * npar,
    BIC = dev.new + log(n) * npar,
    iter = iter
  )
  class(results) = "mrrr"
  return(results)
}
