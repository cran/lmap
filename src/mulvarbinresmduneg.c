//
// Copyright (c) 2020 Frank M.T.A. Busing (e-mail: busing at fsw dot leidenuniv dot nl)
// FreeBSD or 2-Clause BSD or BSD-2 License applies, see Http://www.freebsd.org/copyright/freebsd-license.html
// This is a permissive non-copyleft free software license that is compatible with the GNU GPL. 
//

#include "flib.h"
#include "fmdu.h"
#include "lmap.h"

double mulvarbinresmduneg( const size_t n, const size_t r, double** y, const size_t m, const size_t pu, double** xu, double** bu, const size_t pv, double** xv, double** bv, const bool mains, double* mu, const size_t MAXINNER, const double FCRIT, const size_t MAXITER, const double DCRIT, size_t* lastiter, double* lastdif )
// mulvarbinresmduneg() performs multivariate binary row restricted multidimensional unfolding.
{
  // constants
  const double EPS = DBL_EPSILON;   // 2.2204460492503131e-16
  const double TOL = sqrt( EPS );   // 1.4901161193847656e-08
  const double CRIT = sqrt( TOL );  // 0.00012207031250000000

  // allocate memory
  double** u = getmatrix( n, m, 0.0 );
  double** v = getmatrix( r, m, 0.0 );
  double** q = getmatrix( n, r, 0.0 );
  double** d = getmatrix( n, r, 0.0 );
  double** theta = getmatrix( n, r, 0.0 );
  double** z = getmatrix( n, r, 0.0 );
  double** delta = getmatrix( n, r, 0.0 );

  // compute initial q = 2 y - 1
  for ( size_t i = 1; i <= n; i++ ) {
    for ( size_t j = 1; j <= r; j++ ) q[i][j] = 2.0 * y[i][j] - 1.0;
  }

  // compute initial mu = column means of 4Q
  if ( mains ) {
    for ( size_t j = 1; j <= r; j++ ) {
      double work = 0.0;
      for ( size_t i = 1; i <= n; i++ ) work += q[i][j];
      mu[j] = 4.0 * work / ( double ) ( n );
    }
  }

  // compute initial deviance
  dgemm( false, false, n, m, pu, 1.0, xu, bu, 0.0, u );
  dgemm( false, false, r, m, pv, 1.0, xv, bv, 0.0, v );
  euclidean2( n, m, u, r, v, d );
  for ( size_t j = 1; j <= r; j++ ) {
    double work = mu[j];
    for ( size_t i = 1; i <= n; i++ ) theta[i][j] = work - d[i][j];
  }
  double dold = 0.0;
  for ( size_t i = 1; i <= n; i++ ) {
    for ( size_t j = 1; j <= r; j++ ) {
      const double work = q[i][j] * theta[i][j];
      dold += logl( plogis( work ) );
    }
  }
  dold *= -2.0;

  // start iterations
  size_t iter = 0;
  double ddif = 0.0;
  double dnew = 0.0;
  for ( iter = 1; iter <= MAXITER; iter++ ) {

    // compute working response
    for ( size_t i = 1; i <= n; i++ ) {
      for ( size_t j = 1; j <= r; j++ ) {
        const double work = q[i][j] * theta[i][j];
        z[i][j] = theta[i][j] + 4.0 * q[i][j] * ( 1.0 - plogis( work ) );
      }
    }

    // update main effects
    if ( mains ) {
      for ( size_t j = 1; j <= r; j++ ) {
        double work = 0.0;
        for ( size_t i = 1; i <= n; i++ ) work += z[i][j] + d[i][j];
        mu[j] = work / ( double ) ( n );
      }
    }

    // update U and V
    for ( size_t j = 1; j <= r; j++ ) {
      const double work = mu[j];
      for ( size_t i = 1; i <= n; i++ ) delta[i][j] = work - z[i][j];
    }

    // nonnegative unfolding
    size_t inner = 0;
    double fdif = 0.0;
    resmduneg( n, r, delta, m, pu, xu, bu, pv, xv, bv, d, MAXINNER, FCRIT, &inner, &fdif, false );
    if ( fdif < -1.0 * CRIT ) break;

    // compute deviance
    for ( size_t j = 1; j <= r; j++ ) {
      double work = mu[j];
      for ( size_t i = 1; i <= n; i++ ) theta[i][j] = work - d[i][j];
    }
    dnew = 0.0;
    for ( size_t i = 1; i <= n; i++ ) {
      for ( size_t j = 1; j <= r; j++ ) {
        const double work = q[i][j] * theta[i][j];
        dnew += logl( plogis( work ) );
      }
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
  dgemm( false, false, n, m, pu, 1.0, xu, bu, 0.0, u );
  dgemm( false, false, r, m, pv, 1.0, xv, bv, 0.0, v );
  rotateplusplus( n, m, u, pu, bu, pv, bv );

  // de-allocate memory
  freematrix( u );
  freematrix( v );
  freematrix( q );
  freematrix( d );
  freematrix( theta );
  freematrix( z );
  freematrix( delta );

  return( dnew );
} // mulvarbinresmduneg

void Cmulvarbinresmduneg( int* rn, int* rr, double* ry, int* rm, int* rpu, double* rxu, double* rbu, int* rpv, double* rxv, double* rbv, int* rmains, double* rmu, int* rmaxinner, double* rfcrit, int* rmaxiter, double* rdcrit, double* rdeviance )
// Cmulvarbinresmduneg() performs multivariate binary row restricted multidimensional unfolding.
{
  // transfer to C
  const size_t n = *rn;
  const size_t r = *rr;
  const size_t m = *rm;
  const size_t pu = *rpu;
  const size_t pv = *rpv;
  double** y = getmatrix( n, r, 0.0 );
  for ( size_t j = 1, k = 0; j <= r; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) y[i][j] = ry[k];
  double** xu = getmatrix( n, pu, 0.0 );
  for ( size_t j = 1, k = 0; j <= pu; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) xu[i][j] = rxu[k];
  double** bu = getmatrix( pu, m, 0.0 );
  for ( size_t j = 1, k = 0; j <= m; j++ ) for ( size_t i = 1; i <= pu; i++, k++ ) bu[i][j] = rbu[k];
  double** xv = getmatrix( r, pv, 0.0 );
  for ( size_t j = 1, k = 0; j <= pv; j++ ) for ( size_t i = 1; i <= r; i++, k++ ) xv[i][j] = rxv[k];
  double** bv = getmatrix( pv, m, 0.0 );
  for ( size_t j = 1, k = 0; j <= m; j++ ) for ( size_t i = 1; i <= pv; i++, k++ ) bv[i][j] = rbv[k];
  const bool mains = ( *rmains ) != 0;
  double* mu = getvector( r, 0.0 );
  const size_t MAXINNER = *rmaxinner;
  const double FCRIT = *rfcrit;
  const size_t MAXITER = *rmaxiter;
  const double DCRIT = *rdcrit;

  // analysis
  size_t lastiter = 0;
  double lastdif = 0.0;
  const double dnew = mulvarbinresmduneg( n, r, y, m, pu, xu, bu, pv, xv, bv, mains, mu, MAXINNER, FCRIT, MAXITER, DCRIT, &lastiter, &lastdif );

  // transfer to R
  for ( size_t j = 1, k = 0; j <= m; j++ ) for ( size_t i = 1; i <= pu; i++, k++ ) rbu[k] = bu[i][j];
  for ( size_t j = 1, k = 0; j <= m; j++ ) for ( size_t i = 1; i <= pv; i++, k++ ) rbv[k] = bv[i][j];
  for ( size_t i = 1, k = 0; i <= r; i++, k++ ) rmu[k] = mu[i];
  ( *rmaxiter ) = ( int ) ( lastiter );
  ( *rdcrit ) = lastdif;
  ( *rdeviance ) = dnew;

  // de-allocate memory
  freematrix( y );
  freematrix( xu );
  freematrix( bu );
  freematrix( xv );
  freematrix( bv );
  freevector( mu );

} // Cmulvarbinresmduneg
