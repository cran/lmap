#' Multinomial Restricted MDU
#'
#' The function mru performs multinomial restricted unfolding for a nominal response
#' variable and a set of predictor variables.
#'
#' @param y An N vector of the responses (categorical) or an indicator matrix of size N x C (obligatory when Z is used)
#' @param X An N by P matrix with predictor variables
#' @param Z Design matrix for the class points (V)
#' @param S Positive number indicating the dimensionality of the solution
#' @param start Type of starting values (da: discriminant analysis, random or list with B and V)
#' @param maxiter maximum number of iterations
#' @param dcrit convergence criterion
#' @return Y Matrix Y from input
#' @return Xoriginal Matrix X from input
#' @return X Scaled X matrix
#' @return G class indicator matrix
#' @return ynames class names of response variable
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
#' output = mru(y = y, X = X, S = 2)
#' }
#'
#' @importFrom nnet class.ind
#' @importFrom stats runif
#'
#' @export

mru = function(y, X, Z = NULL, S = 2, start = "da", maxiter = 65536, dcrit = 1e-6)
{

  if(!is.matrix(y)){G = class.ind(y)}
  else if(ncol(y) == 1){G = class.ind(y)}
  else{
    G = y
    y = apply(G, 1, which.max)
  }

  N = nrow(X)
  P = ncol(X)
  C = ncol(G)

  Xoriginal = X
  outx = procx(X)
  X = outx$X
  mx = outx$mx
  sdx = outx$sdx

  # initialization
  if (is.list(start)){
    Bx = start$B
    if(is.null(Z)){V = start$V}
    else{Bz = start$Bz}
  }
  else if( start == "random" ) {
    Bx <- matrix( runif( P * S ), P, S )
    V <- matrix( runif( C * S ), C, S )
  }
  else if ( start == "da" ) {
    tGGinv <- solve( t( G ) %*% G )
    U <- t( X ) %*% G %*% tGGinv
    e <- eigen( ( 1 / N ) * t( X ) %*% X )
    Tmp <- e$vectors %*% diag( sqrt( 1 / e$values ) )
    A <- t( U ) %*% Tmp
    s = svd( A, nu = S, nv = S )
    sV <- matrix( s$v, P, S )
    Bx <- Tmp %*% sV
    V <- tGGinv %*% t( G ) %*% X %*% Bx
  }

  # call C-code
  if(is.null(Z)){
    res = fastmru(G = G, X = X, B = Bx, Z = V, DCRIT = dcrit)
    Bz = NULL
    V = res$V
	npar = S * (P + C - (S-1)/2)
  }
  else{
    Z = cbind(1, Z)
    Bz = solve( t(Z) %*% Z ) %*% t(Z) %*% V
    res = fastmru(G = G, X = X, B = Bx, Z = Z, C = Bz, DCRIT = dcrit)
    Bz = res$C
    V = Z %*% Bz
	npar = S * (P + ncol(Z) - (S-1)/2)
  }

  # make output object
  output = list(
    y = y,
    Xoriginal = Xoriginal,
    G = G,
    X = X,
    ynames = colnames(G),
    xnames = colnames(X),
    mx = mx,
    sdx = sdx,
    U = X %*% res$B,
    Bx = res$B,
    Bz = Bz,
    V = V,
    iter = res$iters,
    deviance = res$deviance,
    npar = npar,
    AIC = res$deviance + 2 * npar,
    BIC = res$deviance + log(N) * npar
)
  class(output) = "mru"
  return(output)
}
