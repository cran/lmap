
#include "library.h"
#include "mdu.h"
#include "lmdu.h"

void Cmulnomrowresmduneg( int* rn, int* rc, double* rg, int* rp, double* rx, int* rm, double* rb, double* rv, double* ru, double* rtheta, int* rmaxiter, double* rdcrit, int* rmaxinner, double* rfcrit, double* rdeviance )
// Function Cmulnomrowresmduneg() performs multinomial row restricted unfolding allowing negative dissimilarities.
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
  size_t p = *rp;
  size_t c = *rc;
  size_t m = *rm;
  double** g = getmatrix( n, c, 0.0 );
  for ( size_t j = 1, k = 0; j <= c; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) g[i][j] = rg[k];
  double** x = getmatrix( n, p, 0.0 );
  for ( size_t j = 1, k = 0; j <= p; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) x[i][j] = rx[k];
  double** b = getmatrix( p, m, 0.0 );
  for ( size_t j = 1, k = 0; j <= m; j++ ) for ( size_t i = 1; i <= p; i++, k++ ) b[i][j] = rb[k];
  double** v = getmatrix( c, m, 0.0 );
  for ( size_t j = 1, k = 0; j <= m; j++ ) for ( size_t i = 1; i <= c; i++, k++ ) v[i][j] = rv[k];
  const size_t MAXITER = *rmaxiter;
  const long double DCRIT = *rdcrit;
  const size_t MAXINNER = *rmaxinner;
  const long double FCRIT = *rfcrit;

  // allocate memory
  double** u = getmatrix( n, m, 0.0 );      // return
  double** theta = getmatrix( n, c, 0.0 );  // return
  double** pi = getmatrix( n, c, 0.0 );     // help
  double** delta = getmatrix( n, c, 0.0 );  // help

  // initialization: u, based on x and b, and pi, based on theta
  gemm( false, false, n, m, p, 1.0L, x, b, 0.0L, u );
  euclidean( n, m, u, c, v, theta );
  for ( size_t i = 1; i <= n; i++ ) {
    long double sum = 0.0L;
    for ( size_t j = 1; j <= c; j++ ) sum += pi[i][j] = exp( -1.0L * theta[i][j] );
    for ( size_t j = 1; j <= c; j++ ) pi[i][j] /= sum;
  }

  // compute old deviance
  long double dold = 0.0L;
  for ( size_t i = 1; i <= n; i++ ) {
    for ( size_t j = 1; j <= c; j++ ) dold += g[i][j] * logl( pi[i][j] );
  }
  dold *= -2.0L;

  // start iterations
  long double ddif = 0.0L;
  long double fdif = 0.0L;
  long double dnew = 0.0L;
  size_t iter = 0;
  size_t inner = 0;
  for ( iter = 1; iter <= MAXITER; iter++ ) {

    // compute delta based on current distances (theta) and pi 
    //bool negs = false;
    for ( size_t i = 1; i <= n; i++ ) {
      for ( size_t j = 1; j <= c; j++ ) {
        delta[i][j] = theta[i][j] - 4.0L * ( g[i][j] - pi[i][j] );
        //if ( delta[i][j] < 0.0 ) negs = true;
      }
    }

    // row restricted mdu allowing negative dissimilarities or not (faster)
    rowresmduneg( n, c, delta, m, p, x, b, v, theta, MAXINNER, FCRIT, &inner, &fdif );
    //if ( negs == true ) rowresmduneg( n, c, delta, m, p, x, b, v, theta, MAXINNER, FCRIT, &inner, &fdif );
    //else rowresmdu( n, c, delta, m, p, x, b, v, theta, MAXINNER, FCRIT, &inner, &fdif );
    if ( fdif < -1.0L * CRIT ) break;

    // compute new pi
    for ( size_t i = 1; i <= n; i++ ) {
      long double sum = 0.0L;
      for ( size_t j = 1; j <= c; j++ ) sum += pi[i][j] = exp( -1.0L * theta[i][j] );
      for ( size_t j = 1; j <= c; j++ ) pi[i][j] /= sum;
    }

    // compute new deviance
    dnew = 0.0L;
    for ( size_t i = 1; i <= n; i++ ) {
      for ( size_t j = 1; j <= c; j++ ) dnew += g[i][j] * logl( pi[i][j] );
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
  gemm( false, false, n, m, p, 1.0, x, b, 0.0, u );
  rotateplusplus( n, m, u, c, v, p, b );

  // transfer to R
  for ( size_t j = 1, k = 0; j <= m; j++ ) for ( size_t i = 1; i <= p; i++, k++ ) rb[k] = b[i][j];
  for ( size_t j = 1, k = 0; j <= m; j++ ) for ( size_t i = 1; i <= c; i++, k++ ) rv[k] = v[i][j];
  for ( size_t j = 1, k = 0; j <= m; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) ru[k] = u[i][j];
  for ( size_t j = 1, k = 0; j <= c; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) rtheta[k] = theta[i][j];
  ( *rmaxiter ) = ( int ) ( iter );
  ( *rdcrit ) = ddif;
  ( *rmaxinner ) = ( int ) ( inner );
  ( *rfcrit ) = fdif;
  ( *rdeviance ) = dnew;

  // de-allocate memory
  freematrix( g );
  freematrix( x );
  freematrix( b );
  freematrix( v );
  freematrix( u );
  freematrix( theta );
  freematrix( pi );
  freematrix( delta );

} // Cmulnomrowresrmduneg
