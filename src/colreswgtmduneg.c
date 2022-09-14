
#include "library.h"
#include "mdu.h"

long double colreswgtmduneg( const size_t n, const size_t m, double** delta, double** w, const size_t p, double** x, const size_t h, double** q, double** b, double** d, const size_t MAXITER, const long double FCRIT, size_t* lastiter, long double* lastdif )
{
  const long double EPS = LDBL_EPSILON;                                              // 2.2204460492503131e-16
  const long double TOL = sqrtl( EPS );                                              // 1.4901161193847656e-08
  const long double CRIT = sqrtl( TOL );                                             // 0.00012207031250000000
  const long double TINY = powl( 10.0L, ( log10l( EPS ) + log10l( TOL ) ) / 2.0L );  // 1.8189894035458617e-12
  const long double DISCRIT = TINY;
  const long double EPSCRIT = 0.25L * TINY;

  // allocate memory
  double** y = getmatrix( m, p, 0.0 );
  double** imb = getmatrix( n, m, 0.0 );
  double** imw = getmatrix( n, m, 0.0 );
  double** xtilde = getmatrix( n, p, 0.0 );
  double** ytilde = getmatrix( m, p, 0.0 );
  double* wc = getvector( m, 0.0 );
  double** hhh = getmatrix( h, h, 0.0 );
  double** hhn = getmatrix( h, n, 0.0 );
  double** hhp = getmatrix( h, p, 0.0 );
  double** hnp = getmatrix( n, p, 0.0 );

  // initialization
  long double scale = 0.0L;
  for ( size_t i = 1; i <= n; i++ ) {
    for ( size_t j = 1; j <= m; j++ ) {
      const long double work = delta[i][j];
      scale += w[i][j] * work * work;
    }
  }

  // update distances and calculate normalized stress
  gemm( false, false, m, p, h, 1.0L, q, b, 0.0L, y );
  euclidean( n, p, x, m, y, d );
  long double fold = 0.0L;
  for ( size_t i = 1; i <= n; i++ ) {
    for ( size_t j = 1; j <= m; j++ ) {
      const long double work = delta[i][j] - d[i][j];
      fold += w[i][j] * work * work;
    }
  }
  fold /= scale;
  long double fnew = 0.0L;

  // start iterations
  size_t iter = 0;
  for ( iter = 1; iter <= MAXITER; iter++ ) {

    // compute original B and W matrices, based on Heiser (1989)
    for ( size_t i = 1; i <= n; i++ ) {
      for ( size_t j = 1; j <= m; j++ ) {
        imb[i][j] = ( delta[i][j] < 0.0 || d[i][j] < DISCRIT ? 0.0 : w[i][j] * delta[i][j] / d[i][j] );
        if ( delta[i][j] < 0.0 ) {
          if ( d[i][j] < DISCRIT ) {
            long double work = fabsl( delta[i][j] );
            imw[i][j] = w[i][j] * ( EPSCRIT + work * work ) / EPSCRIT;
          }
          else imw[i][j] = w[i][j] * ( d[i][j] + fabs( delta[i][j] ) ) / d[i][j];
        }
        else imw[i][j] = w[i][j];
      }
    }
    for ( size_t j = 1; j <= m; j++ ) {
      long double work = 0.0L;
      for ( size_t i = 1; i <= n; i++ ) work += imw[i][j];
      wc[j] = work;
    }

    // compute preliminary updates: xtilde and ytilde
    for ( size_t i = 1; i <= n; i++ ) {
      long double rsb = 0.0L;
      for ( size_t k = 1; k <= m; k++ ) rsb += imb[i][k];
      for ( size_t j = 1; j <= p; j++ ) {
        long double work = 0.0L;
        for ( size_t k = 1; k <= m; k++ ) work += imb[i][k] * y[k][j];
        xtilde[i][j] = rsb * x[i][j] - work;
      }
    }
    for ( size_t i = 1; i <= m; i++ ) {
      long double csb = 0.0L;
      for ( size_t k = 1; k <= n; k++ ) csb += imb[k][i];
      for ( size_t j = 1; j <= p; j++ ) {
        long double work = 0.0L;
        for ( size_t k = 1; k <= n; k++ ) work += imb[k][i] * x[k][j];
        ytilde[i][j] = csb * y[i][j] - work;
      }
    }

    // update x
    gemm( false, false, m, p, n, 1.0L, w, y, 0.0L, hnp );
    for ( size_t i = 1; i <= n; i++ ) {
      long double rsw = 0.0L;
      for ( size_t j = 1; j <= m; j++ ) rsw += imw[i][j];
      for ( size_t j = 1; j <= p; j++ ) x[i][j] = ( xtilde[i][j] + hnp[i][j] ) / rsw;
    }

    // update b
    for ( size_t i = 1; i <= h; i++ ) {
      for ( size_t j = 1; j <= h; j++ ) {
        long double work = 0.0L;
        for ( size_t k = 1; k <= m; k++ ) work += q[k][i] * wc[k] * q[k][j];
        hhh[i][j] = work;
      }
    }
    inverse( h, hhh );
    gemm( true, true, h, n, m, 1.0L, q, imw, 0.0L, hhn );
    gemm( false, false, h, p, n, 1.0L, hhn, x, 0.0L, hhp );
    for ( size_t i = 1; i <= h; i++ ) {
      for ( size_t j = 1; j <= p; j++ ) {
        long double work = 0.0L;
        for ( size_t k = 1; k <= m; k++ ) work += q[k][i] * ytilde[k][j];
        hhp[i][j] += work;
      }
    }
    gemm( false, false, h, p, h, 1.0L, hhh, hhp, 0.0L, b );

    // update y      
    gemm( false, false, m, p, h, 1.0L, q, b, 0.0L, y );

    // update distances and calculate normalized stress
    euclidean( n, p, x, m, y, d );
    fnew = 0.0L;
    for ( size_t i = 1; i <= n; i++ ) {
      for ( size_t j = 1; j <= m; j++ ) {
        long double work = delta[i][j] - d[i][j];
        fnew += w[i][j] * work * work;
      }
    }
    fnew /= scale;

    // check convergence
    ( *lastdif ) = fold - fnew;
    if ( ( *lastdif ) <= -1.0L * CRIT ) break;
    long double fdif = 2.0L * ( *lastdif ) / ( fold + fnew );
    if ( fdif <= FCRIT ) break;
    fold = fnew;
  }
  ( *lastiter ) = iter;

  // rotate to principal axes of x
  rotateplus( n, p, x, h, b );

  // de-allocate memory
  free( ++y[1] ); free( ++y );
  free( ++imb[1] ); free( ++imb );
  free( ++imw[1] ); free( ++imw );
  free( ++xtilde[1] ); free( ++xtilde );
  free( ++ytilde[1] ); free( ++ytilde );
  free( ++wc );
  free( ++hhh[1] ); free( ++hhh );
  free( ++hhn[1] ); free( ++hhn );
  free( ++hhp[1] ); free( ++hhp );
  free( ++hnp[1] ); free( ++hnp );

  return( fnew );
} // colreswgtmduneg

void Ccolreswgtmduneg( int* rn, int* rm, double* rdelta, double* rw, int* rp, double* rx, int* rh, double* rq, double* rb, double* rd, int* rmaxiter, double* rfdif, double* rfvalue )
// Function Ccolresmduneg() performs column restricted multidimensional unfolding allowing negative dissimilarities.
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
  // transfer to C
  size_t n = *rn;
  size_t m = *rm;
  size_t h = *rh;
  size_t p = *rp;
  size_t MAXITER = *rmaxiter;
  double** delta = getmatrix( n, m, 0.0 );
  for ( size_t j = 1, k = 0; j <= m; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) delta[i][j] = rdelta[k];
  double** w = getmatrix( n, m, 0.0 );
  for ( size_t j = 1, k = 0; j <= m; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) w[i][j] = rw[k];
  double** x = getmatrix( n, p, 0.0 );
  for ( size_t j = 1, k = 0; j <= p; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) x[i][j] = rx[k];
  double** q = getmatrix( n, h, 0.0 );
  for ( size_t j = 1, k = 0; j <= h; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) q[i][j] = rq[k];
  double** b = getmatrix( h, p, 0.0 );
  for ( size_t j = 1, k = 0; j <= p; j++ ) for ( size_t i = 1; i <= h; i++, k++ ) b[i][j] = rb[k];
  double** d = getmatrix( n, m, 0.0 );
  long double FCRIT = *rfvalue;

  // run function
  size_t lastiter = 0;
  long double lastdif = 0.0L;
  long double fvalue = colreswgtmduneg( n, m, delta, w, p, x, h, q, b, d, MAXITER, FCRIT, &lastiter, &lastdif );

  // transfer to R
  for ( size_t j = 1, k = 0; j <= p; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) rx[k] = x[i][j];
  for ( size_t j = 1, k = 0; j <= h; j++ ) for ( size_t i = 1; i <= m; i++, k++ ) rq[k] = q[i][j];
  for ( size_t j = 1, k = 0; j <= p; j++ ) for ( size_t i = 1; i <= h; i++, k++ ) rb[k] = b[i][j];
  for ( size_t j = 1, k = 0; j <= m; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) rd[k] = d[i][j];
  ( *rmaxiter ) = ( int ) ( lastiter );
  ( *rfdif ) = lastdif;
  ( *rfvalue ) = fvalue;

  // de-allocate memory
  free( ++delta[1] ); free( ++delta );
  free( ++w[1] ); free( ++w );
  free( ++x[1] ); free( ++x );
  free( ++q[1] ); free( ++q );
  free( ++b[1] ); free( ++b );
  free( ++d[1] ); free( ++d );

} // Ccolreswgtmduneg
