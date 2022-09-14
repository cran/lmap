#' The function lmdu performs logistic mdu with or without predictors to obtain
#' a unsupervised or supervised mapping of binary response variables.
#'
#' This function runs:
#' logistic multidimensional unfolding (if X = NULL)
#' logistic restricted multidimensional unfolding (if X != NULL)
#'
#' @param Y An N times R binary matrix  .
#' @param f Vector with frequencies of response patterns in Y (only applicable if (X = NULL))
#' @param X An N by P matrix with predictor variables
#' @param S Positive number indicating the dimensionality of the solution
#' @param start Either user provided starting values (start should be a list with U and V) or a way to compute starting values (choices: random, svd, ca)
#' @param maxiter maximum number of iterations
#' @param dcrit convergence criterion
#' @return deviance
#'
#' @examples
#' data(dataExample_lmdu)
#' Y = as.matrix(dataExample_lmdu[1:20 , 1:8])
#' X = as.matrix(dataExample_lmdu[1:20 , 9:13])
#' # unsupervised
#' output = lmdu(Y = Y, S = 2)
#'
#' @import tidyverse
#' @import dplyr
#' @importFrom stats plogis rnorm
#'
#' @export
lmdu <- function(Y, f = NULL, X = NULL, S = 2, start = "svd", maxiter = 65536, dcrit = 1e-5)
{


  # checks
  if ( is.null( Y ) ) stop( "missing response variable matrix Y")
  if( ! is.null(X)) if ( nrow(Y) != nrow(X) ) stop( "number of rows in X and Y should match")
  if( S < 1) stop("dimensionality (S) should be larger than 1")
  if( maxiter < 1) stop("maximum number of iterations should be a positive number")
  if( dcrit < 0) stop("converence criterion should be a positive number")
  if( !is.null(f) && !is.null(X) ) stop("f and X cannot both be non-NULL")

  ynames = colnames(Y)

  ###############################################################
  # unsupervised analysis
  ###############################################################
  if( is.null(X) ){
    Xoriginal = NULL
    mx = NULL
    sdx = NULL
    B = NULL
    xnames = NULL

    Yoriginal = Y

    if(is.null(f)){
      R = ncol(Y)
      Ydf = as.data.frame(Y)
      Yr = Ydf %>% group_by_all() %>% summarise(count = n(), .groups = "drop")
      f = Yr$count
      Y = as.matrix(Yr[, 1:R])
      if(sum(Y[1, ]) == 0){ Y = Y[-1, ]; f = f[-1]}
    }

    # Data Coding
    Q = 2 * as.matrix(Y) - 1
    Q[is.na(Q)] <- 0 # Passive treatment of missing data

    N = nrow(Y)

    # Get starting values
    m = colMeans(4 * Q)
    # random starts
    if(start == "random"){
      U = matrix(rnorm(N*S), N, S) #U = matrix(rnorm(N*M), N, S)
      V = matrix(rnorm(R*S), R, S) #V = matrix(rnorm(R*M), R, S)
    }
    else if(start == "svd"){
      Z = as.matrix(outer(rep(1, N), m) + 4 * Q * (1 - plogis(Q * outer(rep(1, N), m))))
      udv = svd(diag(sqrt(f)) %*% scale(Z, center = m, scale = FALSE))
      U = diag(sqrt(1/f)) %*% matrix(udv$u[, 1:S], N, S) %*% diag(udv$d[1:S], nrow = S, ncol = S)
      V = matrix(udv$v[, 1:S], R, S)
    }
    else if(start == "ca"){
      Z = Y/sum(Y)
      r = rowSums(Z)
      c = colSums(Z)
      udv = svd(diag(1/sqrt(r)) %*% (Z - outer(r,c)) %*% diag(1/sqrt(c)))
      U = matrix(udv$u[, 1:S], N, S) %*% diag(udv$d[1:S], nrow = S, ncol = S)
      V = matrix(udv$v[, 1:S], R, S)
    }
    else if( is.list(start) ){
      U = start$U
      if( nrow(U) != N || ncol(U) != S) stop("starting values for U do not have the correct size")
      V = start$V
      if( nrow(V) != R || ncol(V) != S) stop("starting values for V do not have the correct size")
    }
    else stop("no (way of getting) starting values are provided")

    # call C code
    res <- fastmbu( Y, W = f, U, NULL, V, NULL, mains = TRUE, MAXITER = maxiter, DCRIT = dcrit, MAXINNER = 30, FCRIT = 0.000000000001 )
    theta = 0
  }

  ###############################################################
  # supervised analysis
  ###############################################################
  if( !is.null(X) ){

    xnames = colnames(X)
    Yoriginal = Y

    # Data Coding
    Q = 2 * as.matrix(Y) - 1
    Q[is.na(Q)] <- 0 # Passive treatment of missing data

    N = nrow(Q)
    R = ncol(Q)
    P = ncol(X)

    Xoriginal = X
    X = as.matrix(scale(X, center = TRUE, scale = TRUE))
    mx = attr(X, "scaled:center")
    sdx = attr(X, "scaled:scale")

    # starting values m
    m = colMeans(4 * Q)

    # starting values B and V
    if(start == "random"){
      B = matrix(rnorm(P*S), P, S)
      U = X %*% B
      V = matrix(rnorm(R*S), R, S)
    }
    else if(start == "svd"){
      eig.out = eigen(t(X) %*% X)
      iRx = eig.out$vectors %*% diag(1/sqrt(eig.out$values)) %*% t(eig.out$vectors)
      iRxX = iRx %*% t(X)
      Z = as.matrix(outer(rep(1, N), m) + 4 * Q * (1 - plogis(Q * outer(rep(1, N), m))))
      udv = svd(iRxX %*% scale(Z, center = m, scale = FALSE))
      B = iRx %*% matrix(udv$u[, 1:S], P, S) * sqrt(N)
      U = X %*% B
      V = matrix(udv$v[, 1:S], R, S) %*% diag(udv$d[1:S], nrow = S, ncol = S) / sqrt(N)
    }
    else if( is.list(start)){
      B = start$B
      if( nrow(B) != P || ncol(B) != S) stop("starting values for B do not have the correct size")
      V = start$V
      if( nrow(V) != R || ncol(V) != S) stop("starting values for V do not have the correct size")
    }

    # call C-code
    res <- fastmbu( Y, W = NULL, X, B, V, NULL, mains = TRUE, MAXITER = maxiter, DCRIT = dcrit, MAXINNER = 30, FCRIT = 0.000000000001 )
    theta = 0
  }

  # create output object of class "lmdu
  output = list(
    Yoriginal = Yoriginal,
    Xoriginal = Xoriginal,
    Y = Y,
    f = f,
    X = X,
    ynames = ynames,
    xnames = xnames,
    mx = mx,
    sdx = sdx,
    P = plogis(theta),
    m = res$mu,
    U = res$U,
    B = res$BU,
    V = res$V,
    iter = res$iters,
    deviance = res$deviance
  )
  class(output) = "lmdu"
  return(output)
}
