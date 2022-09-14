
#include "library.h"
#include "mdu.h"

long double mduneg( const size_t n, const size_t m, double** delta, const size_t p, double** x, double** y, double** d, const size_t MAXITER, const long double FCRIT, size_t* lastiter, long double* lastdif )
{
  const long double EPS = LDBL_EPSILON;                                              // 2.2204460492503131e-16
  const long double TOL = sqrtl( EPS );                                              // 1.4901161193847656e-08
  const long double CRIT = sqrtl( TOL );                                             // 0.00012207031250000000
  const long double TINY = powl( 10.0L, ( log10l( EPS ) + log10l( TOL ) ) / 2.0L );  // 1.8189894035458617e-12
  const long double DISCRIT = TINY;
  const long double EPSCRIT = 0.25L * TINY;

  // allocate memory
  double** imb = getmatrix( n, m, 0.0 );
  double** imw = getmatrix( n, m, 0.0 );
  double* wr = getvector( n, 0.0 );
  double* wc = getvector( m, 0.0 );
  double** xtilde = getmatrix( n, p, 0.0 );
  double** ytilde = getmatrix( m, p, 0.0 );

  // initialization
  long double scale = 0.0L;
  for ( size_t i = 1; i <= n; i++ ) {
    for ( size_t j = 1; j <= m; j++ ) {
      const long double work = delta[i][j];
      scale += work * work;
    }
  }

  // update distances and calculate normalized stress
  euclidean( n, p, x, m, y, d );
  long double fold = 0.0L;
  for ( size_t i = 1; i <= n; i++ ) {
    for ( size_t j = 1; j <= m; j++ ) {
      long double work = delta[i][j] - d[i][j];
      fold += work * work;
    }
  }
  fold /= scale;
  long double fnew = 0.0L;

  // start unfolding loop
  size_t iter = 0;
  for ( iter = 1; iter <= MAXITER; iter++ ) {

    // compute original B and W matrices, based on Heiser (1989)
    for ( size_t i = 1; i <= n; i++ ) {
      for ( size_t j = 1; j <= m; j++ ) {
        imb[i][j] = ( delta[i][j] < 0.0 || d[i][j] < DISCRIT ? 0.0 : delta[i][j] / d[i][j] );
        if ( delta[i][j] < 0.0 ) {
          if ( d[i][j] < DISCRIT ) {
            const long double work = fabsl( delta[i][j] );
            imw[i][j] = ( EPSCRIT + work * work ) / EPSCRIT;
          }
          else imw[i][j] = ( d[i][j] + fabs( delta[i][j] ) ) / d[i][j];
        }
        else imw[i][j] = 1.0;
      }
    }
    for ( size_t i = 1; i <= n; i++ ) {
      long double work = 0.0L;
      for ( size_t j = 1; j <= m; j++ ) work += imw[i][j];
      wr[i] = work;
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

    // configuration update: x and y
    for ( size_t i = 1; i <= n; i++ ) {
      for ( size_t k = 1; k <= p; k++ ) {
        long double upper = xtilde[i][k];
        for ( size_t j = 1; j <= m; j++ ) upper += imw[i][j] * y[j][k];
        const long double lower = wr[i];
        if ( isnotzero( lower ) ) x[i][k] = upper / lower;
      }
    }
    for ( size_t j = 1; j <= m; j++ ) {
      for ( size_t k = 1; k <= p; k++ ) {
        long double upper = ytilde[j][k];
        for ( size_t i = 1; i <= n; i++ ) upper += imw[i][j] * x[i][k];
        const long double lower = wc[j];
        if ( isnotzero( lower ) ) y[j][k] = upper / lower;
      }
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
  rotateplus( n, p, x, m, y );

  // de-allocate memory
  free( ++imb[1] ); free( ++imb );
  free( ++imw[1] ); free( ++imw );
  free( ++wr );
  free( ++wc );
  free( ++xtilde[1] ); free( ++xtilde );
  free( ++ytilde[1] ); free( ++ytilde );

  return( fnew );
} // mduneg

void Cmduneg( int* rn, int* rm, double* rdelta, int* rp, double* rx, double* ry, double* rd, int* rmaxiter, double* rfdif, double* rfvalue )
// Function Cmduneg() performs multidimensional unfolding allowing negative dissimilarities.
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
  size_t p = *rp;
  size_t MAXITER = *rmaxiter;
  double** delta = getmatrix( n, m, 0.0 );
  for ( size_t j = 1, k = 0; j <= m; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) delta[i][j] = rdelta[k];
  double** x = getmatrix( n, p, 0.0 );
  for ( size_t j = 1, k = 0; j <= p; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) x[i][j] = rx[k];
  double** y = getmatrix( m, p, 0.0 );
  for ( size_t j = 1, k = 0; j <= p; j++ ) for ( size_t i = 1; i <= m; i++, k++ ) y[i][j] = ry[k];
  double** d = getmatrix( n, m, 0.0 );
  long double FCRIT = *rfvalue;

  // run function
  size_t lastiter = 0;
  long double lastdif = 0.0L;
  long double fvalue = mduneg( n, m, delta, p, x, y, d, MAXITER, FCRIT, &lastiter, &lastdif );

  // transfer to R
  for ( size_t j = 1, k = 0; j <= p; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) rx[k] = x[i][j];
  for ( size_t j = 1, k = 0; j <= p; j++ ) for ( size_t i = 1; i <= m; i++, k++ ) ry[k] = y[i][j];
  for ( size_t j = 1, k = 0; j <= m; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) rd[k] = d[i][j];
  ( *rmaxiter ) = ( int ) ( lastiter );
  ( *rfdif ) = lastdif;
  ( *rfvalue ) = fvalue;

  // de-allocate memory
  free( ++delta[1] ); free( ++delta );
  free( ++x[1] ); free( ++x );
  free( ++y[1] ); free( ++y );
  free( ++d[1] ); free( ++d );

} // Cmduneg
