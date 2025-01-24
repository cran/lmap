#'
#' Fast version of mru. It runs mru without input checks.
#'
#' @param G indicator matrix of the response variable
#' @param X matrix with predictor variables
#' @param B starting values of the regression weights
#' @param Z starting values for class locations
#' @param C matrix with coefficients for class points, V = ZC
#' @param MAXINNER maximum number of iterations in the inner loop
#' @param FCRIT convergence criterion for STRESS in the inner loop
#' @param MAXITER maximum number of iterations in the outer loop
#' @param DCRIT convergence criterion for the deviance
#' @param error.check extensive check validity input parameters (default = FALSE).
#'
#' @return B estimated regression weights
#' @return V estimated class locations
#' @return Lastinner number of iterations in the last call to STRESS
#' @return Lastfdif last difference in STRESS values in the inner loop
#' @return lastouter number of iterations in the outer loop
#' @return lastddif last difference in deviances in outer loop
#' @return deviance obtained deviance
#'
#'
#'
#' @export
#' @useDynLib lmap, .registration=TRUE

fastmru <- function( G = NULL, X = NULL, B = NULL, Z = NULL, C = NULL, MAXINNER = 32, FCRIT = 0.001, MAXITER = 65536, DCRIT = 0.000001, error.check = FALSE ){
  # check for input errors
  if ( error.check == TRUE ) {

    # available
    if ( is.null( G ) ) stop( "missing membership data G")
    if ( is.null( X ) ) stop( "missing row data")
    if ( is.null( Z ) ) stop( "missing column data")


    # G
    if ( !is.matrix( G ) ) stop( "G is not a matrix" )
    if ( !is.numeric( G ) ) stop( "G is not numeric" )
    # add G specifics
    n <- nrow( G )
    c <- ncol( G )

    # X
    if ( !is.matrix( X ) ) stop( "X is not a matrix" )
    if ( !is.numeric( X ) ) stop( "X is not numeric" )
    if ( n != nrow( X ) ) stop( "number of rows G and X do not match")
    if ( is.null( B ) ) {
      px <- 0
      m <- ncol( X )
    }
    else {
      if ( !is.matrix( B ) ) stop( "B is not a matrix" )
      if ( !is.numeric( B ) ) stop( "B is not numeric" )
      if ( ncol( X ) != nrow( B ) ) stop( "number of rows B do not match number of columns X")
      px <- ncol( X )
      m <- ncol( B )
    }

    # Z
    if ( !is.matrix( Z ) ) stop( "Z is not a matrix" )
    if ( !is.numeric( Z ) ) stop( "Z is not numeric" )
    if ( c != nrow( Z ) ) stop( "number of columns G and rows Z do not match")
    if ( is.null( C ) ) {
      if ( m != ncol( Z ) ) stop( "number of columns Z and dimensionality do not match")
      pz <- 0
    }
    else {
      if ( m != ncol( C ) ) stop( "number of columns C and dimensionality do not match")
      if ( !is.matrix( C ) ) stop( "C is not a matrix" )
      if ( !is.numeric( C ) ) stop( "C is not numeric" )
      if ( ncol( Z ) != nrow( C ) ) stop( "number of rows V do not match number of columns V")
      pz <- ncol( Z )
    }

    # MAXINNER
    if ( MAXINNER < 0 ) stop( "negative maximum number of inner iterations MAXINNER not allowed" )

    # FCRIT
    if ( FCRIT < 0.0 ) stop( "negative unfolding function convergence criterion not allowed" )

    # MAXITER
    if ( MAXITER < 0 ) stop( "negative maximum number of iterations MAXITER not allowed" )

    # DCRIT
    if ( DCRIT < 0.0 ) stop( "negative likelihood function convergence criterion not allowed" )
  }

  # initialization
  n <- nrow( G )
  c <- ncol( G )
  m <- ifelse( is.null( B ), ncol( X ), ncol( B ) )
  px <- ifelse( is.null( B ), 0, ncol( X ) )
  m <- ifelse( is.null( C ), ncol( Z ), ncol( C ) )
  pz <- ifelse( is.null( C ), 0, ncol( Z ) )
  U <- matrix( 0, n, m )
  V <- matrix( 0, c, m )
  theta <- matrix( 0, n, c )
  deviance <- 0.0

# execution
  if ( is.null( B ) ) {
    # no restrictions
    if ( is.null( C ) ) stop( "invalid model: need at least restriction for one mode" )
    # column restrictions
    else res <- ( .C( "Cmulnomcolresmduneg",n=as.integer(n),c=as.integer(c),G=as.double(G),m=as.integer(m),U=as.double(U),pz=as.integer(pz),Z=as.double(Z),C=as.double(C),theta=as.double(theta),MAXINNER=as.integer(MAXINNER),FCRIT=as.double(FCRIT),MAXITER=as.integer(MAXITER),DCRIT=as.double(DCRIT),deviance=as.double(deviance), PACKAGE = "lmap" ) )
  }
  else {
    # row restrictions
    if ( is.null( C ) ) res <- ( .C( "Cmulnomrowresmduneg",n=as.integer(n),c=as.integer(c),G=as.double(G),m=as.integer(m),px=as.integer(px),X=as.double(X),B=as.double(B),V=as.double(V),theta=as.double(theta),MAXINNER=as.integer(MAXINNER),FCRIT=as.double(FCRIT),MAXITER=as.integer(MAXITER),DCRIT=as.double(DCRIT),deviance=as.double(deviance), PACKAGE = "lmap" ) )
    # row and column restrictions
    else res <- ( .C( "Cmulnomresmduneg",n=as.integer(n),c=as.integer(c),G=as.double(G),m=as.integer(m),px=as.integer(px),X=as.double(X),B=as.double(B),pz=as.integer(pz),Z=as.double(Z),C=as.double(C),theta=as.double(theta),MAXINNER=as.integer(MAXINNER),FCRIT=as.double(FCRIT),MAXITER=as.integer(MAXITER),DCRIT=as.double(DCRIT),deviance=as.double(deviance), PACKAGE = "lmap" ) )
  }


  # finalization
  if ( is.null( B ) ) U <- matrix( res$X, n, m )
  else {
    B <- matrix( res$B, px, m )
    U <- X %*% B
  }
  if ( is.null( C ) ) V <- matrix( res$V, c, m )
  else {
    C <- matrix( res$C, pz, m )
    V <- Z %*% C
  }

  return( list( B=B,
                U=U,
                C=C,
                V=V,
                lastiter=res$MAXITER,
                lastddif=res$DCRIT,
                deviance=res$deviance ) )

} # fastmru
