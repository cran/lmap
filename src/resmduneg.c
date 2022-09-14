
#include "library.h"
#include "mdu.h"

long double resmduneg( const size_t n, const size_t m, double** delta, const size_t p, const size_t hx, double** qx, double** bx, const size_t hy, double** qy, double** by, double** d, const size_t MAXITER, const long double FCRIT, size_t* lastiter, long double* lastdif )
{
  const long double EPS = LDBL_EPSILON;                                              // 2.2204460492503131e-16
  const long double TOL = sqrtl( EPS );                                              // 1.4901161193847656e-08
  const long double CRIT = sqrtl( TOL );                                             // 0.00012207031250000000
  const long double TINY = powl( 10.0L, ( log10l( EPS ) + log10l( TOL ) ) / 2.0L );  // 1.8189894035458617e-12
  const long double DISCRIT = TINY;
  const long double EPSCRIT = 0.25L * TINY;

  // allocate memory
  double** x = getmatrix( n, p, 0.0 );
  double** y = getmatrix( m, p, 0.0 );
  double** imb = getmatrix( n, m, 0.0 );
  double** imw = getmatrix( n, m, 0.0 );
  double** xtilde = getmatrix( n, p, 0.0 );
  double** ytilde = getmatrix( m, p, 0.0 );
  double* wr = getvector( n, 0.0 );
  double* wc = getvector( m, 0.0 );
  double** hxx = getmatrix( hx, hx, 0.0 );
  double** hhm = getmatrix( hx, m, 0.0 );
  double** hhp = getmatrix( hx, p, 0.0 );
  double** hmp = getmatrix( m, p, 0.0 );
  double** hyy = getmatrix( hy, hy, 0.0 );
  double** hhn = getmatrix( hy, m, 0.0 );
  double** hnp = getmatrix( m, p, 0.0 );

  // initialization
  long double scale = 0.0L;
  for ( size_t i = 1; i <= n; i++ ) {
    for ( size_t j = 1; j <= m; j++ ) {
      const long double work = delta[i][j];
      scale += work * work;
    }
  }

  // update distances and calculate normalized stress
  gemm( false, false, n, p, hx, 1.0L, qx, bx, 0.0L, x );
  gemm( false, false, m, p, hy, 1.0L, qy, by, 0.0L, y );
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

    // update bx
    for ( size_t i = 1; i <= hx; i++ ) {
      for ( size_t j = 1; j <= hx; j++ ) {
        long double work = 0.0L;
        for ( size_t k = 1; k <= n; k++ ) work += qx[k][i] * wr[k] * qx[k][j];
        hxx[i][j] = work;
      }
    }
    inverse( hx, hxx );
    gemm( true, false, hx, m, n, 1.0L, qx, imw, 0.0L, hhm );
    gemm( false, false, hx, p, m, 1.0L, hhm, y, 0.0L, hhp );
    for ( size_t i = 1; i <= hx; i++ ) {
      for ( size_t j = 1; j <= p; j++ ) {
        long double work = 0.0L;
        for ( size_t k = 1; k <= n; k++ ) work += qx[k][i] * xtilde[k][j];
        hhp[i][j] += work;
      }
    }
    gemm( false, false, hx, p, hx, 1.0L, hxx, hhp, 0.0L, bx );

    // update x      
    gemm( false, false, n, p, hx, 1.0L, qx, bx, 0.0L, x );

    // update by
    for ( size_t i = 1; i <= hy; i++ ) {
      for ( size_t j = 1; j <= hy; j++ ) {
        long double work = 0.0L;
        for ( size_t k = 1; k <= m; k++ ) work += qy[k][i] * wc[k] * qy[k][j];
        hyy[i][j] = work;
      }
    }
    inverse( hy, hyy );
    gemm( true, true, hy, n, m, 1.0L, qy, imw, 0.0L, hhn );
    gemm( false, false, hy, p, n, 1.0L, hhn, x, 0.0L, hhp );
    for ( size_t i = 1; i <= hy; i++ ) {
      for ( size_t j = 1; j <= p; j++ ) {
        long double work = 0.0L;
        for ( size_t k = 1; k <= m; k++ ) work += qy[k][i] * ytilde[k][j];
        hhp[i][j] += work;
      }
    }
    gemm( false, false, hy, p, hy, 1.0L, hyy, hhp, 0.0L, by );

    // update y      
    gemm( false, false, m, p, hy, 1.0L, qy, by, 0.0L, y );

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
  rotateplusplus( n, p, x, hx, bx, hy, by );

  // de-allocate memory
  free( ++x[1] ); free( ++x );
  free( ++y[1] ); free( ++y );
  free( ++imb[1] ); free( ++imb );
  free( ++imw[1] ); free( ++imw );
  free( ++xtilde[1] ); free( ++xtilde );
  free( ++ytilde[1] ); free( ++ytilde );
  free( ++wr );
  free( ++wc );
  free( ++hxx[1] ); free( ++hxx );
  free( ++hhm[1] ); free( ++hhm );
  free( ++hhp[1] ); free( ++hhp );
  free( ++hmp[1] ); free( ++hmp );
  free( ++hyy[1] ); free( ++hyy );
  free( ++hhn[1] ); free( ++hhn );
  free( ++hnp[1] ); free( ++hnp );

  return( fnew );
} // resmduneg

void Cresmduneg( int* rn, int* rm, double* rdelta, int* rp, int* rhx, double* rqx, double* rbx, int* rhy, double* rqy, double* rby, double* rd, int* rmaxiter, double* rfdif, double* rfvalue )
// Function Cresmduneg() performs row restricted multidimensional unfolding.
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
  size_t hx = *rhx;
  size_t hy = *rhy;
  size_t p = *rp;
  size_t MAXITER = *rmaxiter;
  double** delta = getmatrix( n, m, 0.0 );
  for ( size_t j = 1, k = 0; j <= m; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) delta[i][j] = rdelta[k];
  double** qx = getmatrix( n, hx, 0.0 );
  for ( size_t j = 1, k = 0; j <= hx; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) qx[i][j] = rqx[k];
  double** bx = getmatrix( hx, p, 0.0 );
  for ( size_t j = 1, k = 0; j <= p; j++ ) for ( size_t i = 1; i <= hx; i++, k++ ) bx[i][j] = rbx[k];
  double** qy = getmatrix( n, hy, 0.0 );
  for ( size_t j = 1, k = 0; j <= hy; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) qy[i][j] = rqy[k];
  double** by = getmatrix( hy, p, 0.0 );
  for ( size_t j = 1, k = 0; j <= p; j++ ) for ( size_t i = 1; i <= hy; i++, k++ ) by[i][j] = rby[k];
  double** d = getmatrix( n, m, 0.0 );
  long double FCRIT = *rfvalue;

  // run function
  size_t lastiter = 0;
  long double lastdif = 0.0L;
  long double fvalue = resmduneg( n, m, delta, p, hx, qx, bx, hy, qy, by, d, MAXITER, FCRIT, &lastiter, &lastdif );

  // transfer to R
  for ( size_t j = 1, k = 0; j <= hx; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) rqx[k] = qx[i][j];
  for ( size_t j = 1, k = 0; j <= p; j++ ) for ( size_t i = 1; i <= hx; i++, k++ ) rbx[k] = bx[i][j];
  for ( size_t j = 1, k = 0; j <= hy; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) rqy[k] = qy[i][j];
  for ( size_t j = 1, k = 0; j <= p; j++ ) for ( size_t i = 1; i <= hy; i++, k++ ) rby[k] = by[i][j];
  for ( size_t j = 1, k = 0; j <= m; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) rd[k] = d[i][j];
  ( *rmaxiter ) = ( int ) ( lastiter );
  ( *rfdif ) = ( int ) ( lastdif );
  ( *rfvalue ) = fvalue;

  // de-allocate memory
  free( ++delta[1] ); free( ++delta );
  free( ++qx[1] ); free( ++qx );
  free( ++bx[1] ); free( ++bx );
  free( ++qy[1] ); free( ++qy );
  free( ++by[1] ); free( ++by );
  free( ++d[1] ); free( ++d );

} // Cresmduneg
