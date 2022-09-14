#'
#' Fast version of mru. It runs mru without input checks.
#'
#' @param G indicator matrix of the response variable
#' @param X matrix with predictor variables
#' @param B starting values of the regression weights
#' @param V starting values for class locations
#' @param MAXITER maximum number of iterations in the outer loop
#' @param DCRIT convergence criterion for the deviance
#' @param MAXINNER maximum number of iterations in the inner loop
#' @param FCRIT convergence criterion for STRESS in the inner loop
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
#' @useDynLib lmap

fastmru <- function( G = NULL, X = NULL, B = NULL, V = NULL, MAXITER = 65536, DCRIT = 0.000001, MAXINNER = 32, FCRIT = 0.001 )
{
  # initialization
  n <- nrow( G )
  c <- ncol( G )
  p <- ncol( X )
  m <- ncol( B )
  U <- matrix( 0, n, m )
  theta <- matrix( 0, n, c )
  deviance <- 0.0

  # execution
  res <- ( .C( "Cmulnomrowresmduneg",
             n=as.integer(n),
             c=as.integer(c),
             G=as.double(G),
             p=as.integer(p),
             X=as.double(X),
             m=as.integer(m),
             B=as.double(B),
             V=as.double(V),
             U=as.double(U),
             theta=as.double(theta),
             MAXITER=as.integer(MAXITER),
             DCRIT=as.double(DCRIT),
             MAXINNER=as.integer(MAXINNER),
             FCRIT=as.double(FCRIT),
             deviance=as.double(deviance), PACKAGE = "lmap" ) )

  # finalization
  B <- matrix( res$B, p, m )
  V <- matrix( res$V, c, m )

  return( list( B=B, V=V,
                lastinner=res$MAXINNER, lastfdif=res$FCRIT,
                lastouter=res$MAXITER, lastddif=res$DCRIT, deviance=res$deviance ) )

} # fastmru
