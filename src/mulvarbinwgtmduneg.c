
#include "library.h"
#include "mdu.h"
#include "lmdu.h"

void Cmulvarbinwgtmduneg( int* rn, int* rr, double* ry, double* rw, int* rm, double* ru, double* rv, int* rmains, double* rmu, int* rmaxiter, double* rdcrit, int* rmaxinner, double* rfcrit, double* rdeviance )
// Function Cmulvarbinwgtmduneg() performs multivariate binary weighted multidimensional unfolding.
// Copyright (C) 2020 Frank M.T.A. Busing (e-mail: busing at fsw dot leidenuniv dot nl)
// This function is free software:
// you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation.
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY;
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
// You should have received a copy of the GNU General Public License along with this function.
// If not, see <https://www.gnu.org/licenses/>.
{
  // constants
  const long double EPS = LDBL_EPSILON;   // 2.2204460492503131e-16
  const long double TOL = sqrtl( EPS );   // 1.4901161193847656e-08
  const long double CRIT = sqrtl( TOL );  // 0.00012207031250000000

  // transfer to C
  size_t n = *rn;
  size_t r = *rr;
  size_t m = *rm;
  bool mains = ( *rmains ) != 0;
  double** y = getmatrix( n, r, 0.0 );
  for ( size_t j = 1, k = 0; j <= r; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) y[i][j] = ry[k];
  double** w = getmatrix( n, r, 0.0 );
  for ( size_t j = 1, k = 0; j <= r; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) w[i][j] = rw[k];
  double** u = getmatrix( n, m, 0.0 );
  for ( size_t j = 1, k = 0; j <= m; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) u[i][j] = ru[k];
  double** v = getmatrix( r, m, 0.0 );
  for ( size_t j = 1, k = 0; j <= m; j++ ) for ( size_t i = 1; i <= r; i++, k++ ) v[i][j] = rv[k];
  size_t MAXITER = *rmaxiter;
  long double DCRIT = *rdcrit;
  size_t MAXINNER = *rmaxinner;
  long double FCRIT = *rfcrit;

  // allocate memory
  double* mu = getvector( r, 0.0 );         // return
  double** q = getmatrix( n, r, 0.0 );      // help
  double** d = getmatrix( n, r, 0.0 );      // help
  double** theta = getmatrix( n, r, 0.0 );  // help
  double** z = getmatrix( n, r, 0.0 );      // help
  double** delta = getmatrix( n, r, 0.0 );  // help
  double* wr = getvector( r, 0.0 );         // help

  // compute initial q = 2 y - 1
  for ( size_t i = 1; i <= n; i++ ) {
    for ( size_t j = 1; j <= r; j++ ) q[i][j] = 2.0 * y[i][j] - 1.0;
  }

  // compute initial mu = column means of 4Q
  if ( mains ) {
    for ( size_t j = 1; j <= r; j++ ) {
      long double work = 0.0L;
      for ( size_t i = 1; i <= n; i++ ) {
        wr[j] += w[i][j];
        work += w[i][j] * q[i][j];
      }
      mu[j] = 4.0L * work / wr[j];
    }
  }

  // compute initial deviance
  euclidean( n, m, u, r, v, d );
  if ( mains == true ) {
    for ( size_t j = 1; j <= r; j++ ) {
      long double work = mu[j];
      for ( size_t i = 1; i <= n; i++ ) theta[i][j] = work - d[i][j];
    }
  }
  else {
    for ( size_t j = 1; j <= r; j++ ) {
      for ( size_t i = 1; i <= n; i++ ) theta[i][j] = -1.0 * d[i][j];
    }
  }
  long double dold = 0.0L;
  for ( size_t i = 1; i <= n; i++ ) {
    for ( size_t j = 1; j <= r; j++ ) {
      const long double work = q[i][j] * theta[i][j];
      dold += w[i][j] * logl( plogis( work ) );
    }
  }
  dold *= -2.0L;
  long double dnew = 0.0L;

  // start iterations
  size_t iter = 0;
  long double ddif = 0.0L;
  size_t inner = 0;
  long double fdif = 0.0L;
  for ( iter = 1; iter <= MAXITER; iter++ ) {

    // compute working response
    for ( size_t i = 1; i <= n; i++ ) {
      for ( size_t j = 1; j <= r; j++ ) {
        const long double work = q[i][j] * theta[i][j];
        z[i][j] = theta[i][j] + 4.0L * q[i][j] * ( 1.0L - plogis( work ) );
      }
    }

    // update main effects
    if ( mains ) {
      for ( size_t j = 1; j <= r; j++ ) {
        long double work = 0.0L;
        for ( size_t i = 1; i <= n; i++ ) work += w[i][j] *( z[i][j] + d[i][j] );
        mu[j] = work / wr[j];
      }

      // update U and V
      for ( size_t j = 1; j <= r; j++ ) {
        const long double work = mu[j];
        for ( size_t i = 1; i <= n; i++ ) delta[i][j] = work - z[i][j];
      }
    }
    else {

      // update U and V
      for ( size_t j = 1; j <= r; j++ ) {
        for ( size_t i = 1; i <= n; i++ ) delta[i][j] = -1.0 * z[i][j];
      }
    }

    // weighted nonnegative unfolding
    wgtmduneg( n, r, delta, w, m, u, v, d, MAXINNER, FCRIT, &inner, &fdif );
    if ( fdif < -1.0L * CRIT ) break;

    // compute deviance
    if ( mains == true ) {
      for ( size_t j = 1; j <= r; j++ ) {
        long double work = mu[j];
        for ( size_t i = 1; i <= n; i++ ) theta[i][j] = work - d[i][j];
      }
    }
    else {
      for ( size_t j = 1; j <= r; j++ ) {
        for ( size_t i = 1; i <= n; i++ ) theta[i][j] = -1.0 * d[i][j];
      }
    }
    dnew = 0.0L;
    for ( size_t i = 1; i <= n; i++ ) {
      for ( size_t j = 1; j <= r; j++ ) {
        const long double work = q[i][j] * theta[i][j];
        dnew += w[i][j] * logl( plogis( work ) );
      }
    }
    dnew *= -2.0L;

    // check convergence
    long double lastdif = dold - dnew;
    if ( lastdif <= -1.0L * CRIT ) break;
    ddif = 2.0L * lastdif / ( dold + dnew );
    if ( ddif <= DCRIT ) break;
    dold = dnew;
  }

  // rotate solution to principal axes
  double* wn = getvector( n, 0.0 );
  for ( size_t i = 1; i <= n; i++ ) {
    long double work = 0.0L;
    for ( size_t j = 1; j <= r; j++ ) work += w[i][j];
    wn[i] = work / ( long double ) ( r );
  }
  weightedrotateplus( n, m, u, wn, r, v );

  // transfer to R
  for ( size_t j = 1, k = 0; j <= m; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) ru[k] = u[i][j];
  for ( size_t j = 1, k = 0; j <= m; j++ ) for ( size_t i = 1; i <= r; i++, k++ ) rv[k] = v[i][j];
  for ( size_t i = 1, k = 0; i <= r; i++, k++ ) rmu[k] = mu[i];
  ( *rmaxiter ) = ( int ) ( iter );
  ( *rdcrit ) = ddif;
  ( *rmaxinner ) = ( int ) ( inner );
  ( *rfcrit ) = fdif;
  ( *rdeviance ) = dnew;

  // de-allocate memory
  freematrix( y );
  freematrix( w );
  freematrix( u );
  freematrix( v );
  freematrix( q );
  freevector( mu );
  freematrix( d );
  freematrix( theta );
  freematrix( z );
  freematrix( delta );
  freevector( wr );
  freevector( wn );

} // Cmulvarbinwgtmduneg
