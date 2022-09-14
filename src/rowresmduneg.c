
#include "library.h"
#include "mdu.h"

long double rowresmduneg( const size_t n, const size_t m, double** delta, const size_t p, const size_t h, double** q, double** b, double** y, double** d, const size_t MAXITER, const long double FCRIT, size_t* lastiter, long double* lastdif )
{
  const long double EPS = LDBL_EPSILON;                                              // 2.2204460492503131e-16
  const long double TOL = sqrtl( EPS );                                              // 1.4901161193847656e-08
  const long double CRIT = sqrtl( TOL );                                             // 0.00012207031250000000
  const long double TINY = powl( 10.0L, ( log10l( EPS ) + log10l( TOL ) ) / 2.0L );  // 1.8189894035458617e-12
  const long double DISCRIT = TINY;
  const long double EPSCRIT = 0.25L * TINY;

  // allocate memory
  double** x = getmatrix( n, p, 0.0 );
  double** imb = getmatrix( n, m, 0.0 );
  double** imw = getmatrix( n, m, 0.0 );
  double** xtilde = getmatrix( n, p, 0.0 );
  double** ytilde = getmatrix( m, p, 0.0 );
  double* wr = getvector( n, 0.0 );
  double** hhh = getmatrix( h, h, 0.0 );
  double** hhm = getmatrix( h, m, 0.0 );
  double** hhp = getmatrix( h, p, 0.0 );
  double** hmp = getmatrix( m, p, 0.0 );

  // initialization
  long double scale = 0.0L;
  for ( size_t i = 1; i <= n; i++ ) {
    for ( size_t j = 1; j <= m; j++ ) {
      const long double work = delta[i][j];
      scale += work * work;
    }
  }

  // update distances and calculate normalized stress
  gemm( false, false, n, p, h, 1.0L, q, b, 0.0L, x );
  euclidean( n, p, x, m, y, d );
  long double fold = 0.0L;
  for ( size_t i = 1; i <= n; i++ ) {
    for ( size_t j = 1; j <= m; j++ ) {
      const long double work = delta[i][j] - d[i][j];
      fold += work * work;
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
        imb[i][j] = ( delta[i][j] < 0.0 || d[i][j] < DISCRIT ? 0.0 : delta[i][j] / d[i][j] );
        if ( delta[i][j] < 0.0 ) {
          if ( d[i][j] < DISCRIT ) {
            long double work = fabsl( delta[i][j] );
            imw[i][j] = ( EPSCRIT + work * work ) / EPSCRIT;
          }
          else imw[i][j] = ( d[i][j] + fabs( delta[i][j] ) ) / d[i][j];
        }
        else imw[i][j] = 1.0;
      }
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

    // update b
    for ( size_t i = 1; i <= n; i++ ) {
      long double work = 0.0L;
      for ( size_t j = 1; j <= m; j++ ) work += imw[i][j];
      wr[i] = work;
    }
    for ( size_t i = 1; i <= h; i++ ) {
      for ( size_t j = 1; j <= h; j++ ) {
        long double work = 0.0L;
        for ( size_t k = 1; k <= n; k++ ) work += q[k][i] * wr[k] * q[k][j];
        hhh[i][j] = work;
      }
    }
    inverse( h, hhh );
    gemm( true, false, h, m, n, 1.0L, q, imw, 0.0L, hhm );
    gemm( false, false, h, p, m, 1.0L, hhm, y, 0.0L, hhp );
    for ( size_t i = 1; i <= h; i++ ) {
      for ( size_t j = 1; j <= p; j++ ) {
        long double work = 0.0L;
        for ( size_t k = 1; k <= n; k++ ) work += q[k][i] * xtilde[k][j];
        hhp[i][j] += work;
      }
    }
    gemm( false, false, h, p, h, 1.0L, hhh, hhp, 0.0L, b );

    // update x      
    gemm( false, false, n, p, h, 1.0L, q, b, 0.0L, x );

    // update y
    gemm( true, false, m, p, n, 1.0L, imw, x, 0.0L, hmp );
    for ( size_t i = 1; i <= m; i++ ) {
      long double csw = 0.0L;
      for ( size_t j = 1; j <= n; j++ ) csw += imw[j][i];
      for ( size_t j = 1; j <= p; j++ ) y[i][j] = ( ytilde[i][j] + hmp[i][j] ) / csw;
    }

    // update distances and calculate normalized stress
    euclidean( n, p, x, m, y, d );
    fnew = 0.0L;
    for ( size_t i = 1; i <= n; i++ ) {
      for ( size_t j = 1; j <= m; j++ ) {
        long double work = delta[i][j] - d[i][j];
        fnew += work * work;
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
  rotateplusplus( n, p, x, h, b, m, y );

  // de-allocate memory
  free( ++x[1] ); free( ++x );
  free( ++imb[1] ); free( ++imb );
  free( ++imw[1] ); free( ++imw );
  free( ++xtilde[1] ); free( ++xtilde );
  free( ++ytilde[1] ); free( ++ytilde );
  free( ++wr );
  free( ++hhh[1] ); free( ++hhh );
  free( ++hhm[1] ); free( ++hhm );
  free( ++hhp[1] ); free( ++hhp );
  free( ++hmp[1] ); free( ++hmp );

  return( fnew );
} // rowresmduneg

void Crowresmduneg( int* rn, int* rm, double* rdelta, int* rp, int* rh, double* rq, double* rb, double* ry, double* rd, int* rmaxiter, double* rfdif, double* rfvalue )
// Function Crowresmduneg() performs restricted multidimensional unfolding allowing negative dissimilarities.
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
  double** q = getmatrix( n, h, 0.0 );
  for ( size_t j = 1, k = 0; j <= h; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) q[i][j] = rq[k];
  double** b = getmatrix( h, p, 0.0 );
  for ( size_t j = 1, k = 0; j <= p; j++ ) for ( size_t i = 1; i <= h; i++, k++ ) b[i][j] = rb[k];
  double** y = getmatrix( m, p, 0.0 );
  for ( size_t j = 1, k = 0; j <= p; j++ ) for ( size_t i = 1; i <= m; i++, k++ ) y[i][j] = ry[k];
  double** d = getmatrix( n, m, 0.0 );
  long double FCRIT = *rfvalue;

  // run function
  size_t lastiter = 0;
  long double lastdif = 0.0L;
  long double fvalue = rowresmduneg( n, m, delta, p, h, q, b, y, d, MAXITER, FCRIT, &lastiter, &lastdif );

  // transfer to R
  for ( size_t j = 1, k = 0; j <= h; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) rq[k] = q[i][j];
  for ( size_t j = 1, k = 0; j <= p; j++ ) for ( size_t i = 1; i <= h; i++, k++ ) rb[k] = b[i][j];
  for ( size_t j = 1, k = 0; j <= p; j++ ) for ( size_t i = 1; i <= m; i++, k++ ) ry[k] = y[i][j];
  for ( size_t j = 1, k = 0; j <= m; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) rd[k] = d[i][j];
  ( *rmaxiter ) = ( int ) ( lastiter );
  ( *rfdif ) = ( int ) ( lastdif );
  ( *rfvalue ) = fvalue;

  // de-allocate memory
  free( ++delta[1] ); free( ++delta );
  free( ++q[1] ); free( ++q );
  free( ++b[1] ); free( ++b );
  free( ++y[1] ); free( ++y );
  free( ++d[1] ); free( ++d );

} // Crowresmduneg
