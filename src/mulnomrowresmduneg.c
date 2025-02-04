//
// Copyright (c) 2020 Frank M.T.A. Busing (e-mail: busing at fsw dot leidenuniv dot nl)
// FreeBSD or 2-Clause BSD or BSD-2 License applies, see Http://www.freebsd.org/copyright/freebsd-license.html
// This is a permissive non-copyleft free software license that is compatible with the GNU GPL. 
//

#include "flib.h"
#include "fmdu.h"
#include "lmap.h"

double mulnomrowresmduneg( const size_t n, const size_t nc, double** g, const size_t m, const size_t px, double** x, double** b, double** v, double** theta, const size_t MAXINNER, const double FCRIT, const size_t MAXITER, const double DCRIT, size_t* lastiter, double* lastdif )
// mulnomrowresmduneg() performs multinomial row restricted unfolding allowing negative dissimilarities.
{
  // constants
  const double EPS = DBL_EPSILON;   // 2.2204460492503131e-16
  const double TOL = sqrt( EPS );   // 1.4901161193847656e-08
  const double CRIT = sqrt( TOL );  // 0.00012207031250000000

  // allocate memory
  double** u = getmatrix( n, m, 0.0 );
  double** pi = getmatrix( n, nc, 0.0 );
  double** delta = getmatrix( n, nc, 0.0 );
  int** fv = getimatrix( nc, m, 0 );

  // initialization: u, based on x and b, and pi, based on theta
  dgemm( false, false, n, m, px, 1.0, x, b, 0.0, u );
  euclidean2( n, m, u, nc, v, theta );
  for ( size_t i = 1; i <= n; i++ ) {
    double sum = 0.0;
    for ( size_t j = 1; j <= nc; j++ ) sum += pi[i][j] = exp( -1.0 * theta[i][j] );
    for ( size_t j = 1; j <= nc; j++ ) pi[i][j] /= sum;
  }

  // compute old deviance
  double dold = 0.0;
  for ( size_t i = 1; i <= n; i++ ) {
    for ( size_t j = 1; j <= nc; j++ ) dold += g[i][j] * logl( pi[i][j] );
  }
  dold *= -2.0;

  // start iterations
  double ddif = 0.0;
  double dnew = 0.0;
  size_t iter = 0;
  for ( iter = 1; iter <= MAXITER; iter++ ) {

    // compute delta based on current distances (theta) and pi 
    //bool negs = false;
    for ( size_t i = 1; i <= n; i++ ) {
      for ( size_t j = 1; j <= nc; j++ ) {
        delta[i][j] = theta[i][j] - 4.0 * ( g[i][j] - pi[i][j] );
        //if ( delta[i][j] < 0.0 ) negs = true;
      }
    }

    // row restricted mdu allowing negative dissimilarities or not (faster)
    size_t inner = 0;
    double fdif = 0.0;
    rowresmduneg( n, nc, delta, m, px, x, b, v, fv, theta, MAXINNER, FCRIT, &inner, &fdif, false );
																											
																						  
    if ( fdif < -1.0 * CRIT ) break;

    // compute new pi
    for ( size_t i = 1; i <= n; i++ ) {
      double sum = 0.0;
      for ( size_t j = 1; j <= nc; j++ ) sum += pi[i][j] = exp( -1.0 * theta[i][j] );
      for ( size_t j = 1; j <= nc; j++ ) pi[i][j] /= sum;
    }

    // compute new deviance
    dnew = 0.0;
    for ( size_t i = 1; i <= n; i++ ) {
      for ( size_t j = 1; j <= nc; j++ ) dnew += g[i][j] * logl( pi[i][j] );
    }
    dnew *= -2.0;

    // check convergence
    ( *lastdif ) = dold - dnew;
    if ( ( *lastdif ) <= -1.0 * CRIT ) break;
    ddif = 2.0 * ( *lastdif ) / ( dold + dnew );
    if ( ddif <= DCRIT ) break;
    dold = dnew;
  }
  ( *lastiter ) = iter;

  // rotate solution to principal axes
  dgemm( false, false, n, m, px, 1.0, x, b, 0.0, u );
  rotateplusplus( n, m, u, nc, v, px, b );

  // de-allocate memory
  freematrix( u );
  freematrix( pi );
  freematrix( delta );
  freeimatrix( fv );

  return( dnew );
} // mulnomrowresrmduneg

void Cmulnomrowresmduneg( int* rn, int* rnc, double* rg, int* rm, int* rpx, double* rx, double* rb, double* rv, double* rtheta, int* rmaxinner, double* rfcrit, int* rmaxiter, double* rdcrit, double* rdeviance )
// Cmulnomrowresmduneg() performs multinomial row restricted unfolding allowing negative dissimilarities.
{
  // transfer to C
  const size_t n = *rn;
  const size_t nc = *rnc;
					   
  const size_t m = *rm;
  const size_t px = *rpx;
  double** g = getmatrix( n, nc, 0.0 );
  for ( size_t j = 1, k = 0; j <= nc; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) g[i][j] = rg[k];
  double** x = getmatrix( n, px, 0.0 );
  for ( size_t j = 1, k = 0; j <= px; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) x[i][j] = rx[k];
  double** b = getmatrix( px, m, 0.0 );
  for ( size_t j = 1, k = 0; j <= m; j++ ) for ( size_t i = 1; i <= px; i++, k++ ) b[i][j] = rb[k];
  double** v = getmatrix( nc, m, 0.0 );
  for ( size_t j = 1, k = 0; j <= m; j++ ) for ( size_t i = 1; i <= nc; i++, k++ ) v[i][j] = rv[k];
  double** u = getmatrix( n, m, 0.0 );
  double** theta = getmatrix( n, nc, 0.0 );
  const size_t MAXITER = *rmaxiter;
  const double DCRIT = *rdcrit;
  const size_t MAXINNER = *rmaxinner;
  const double FCRIT = *rfcrit;

  // analysis
  size_t lastiter = 0;
  double lastdif = 0.0;
  const double dnew = mulnomrowresmduneg( n, nc, g, m, px, x, b, v, theta, MAXINNER, FCRIT, MAXITER, DCRIT, &lastiter, &lastdif );

  // transfer to R
  for ( size_t j = 1, k = 0; j <= m; j++ ) for ( size_t i = 1; i <= px; i++, k++ ) rb[k] = b[i][j];
  for ( size_t j = 1, k = 0; j <= m; j++ ) for ( size_t i = 1; i <= nc; i++, k++ ) rv[k] = v[i][j];
																								  
  for ( size_t j = 1, k = 0; j <= nc; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) rtheta[k] = theta[i][j];
  ( *rmaxiter ) = ( int ) ( lastiter );
  ( *rdcrit ) = lastdif;
  ( *rdeviance ) = dnew;

  // de-allocate memory
  freematrix( g );
  freematrix( x );
  freematrix( b );
  freematrix( v );
  freematrix( u );
  freematrix( theta );

} // Cmulnomrowresrmduneg
