
#include "library.h"

bool isequal( const long double d1, const long double d2 )
// eps hold the minimum distance between the values 
// that will be considered as the numbers are equal
// considering the magnitude of the numbers
{
  const long double EPS = LDBL_EPSILON;
  long double eps1 = fabsl( d1 );
  long double eps2 = fabsl( d2 );
  long double eps = ( eps1 > eps2 ) ? eps1 : eps2;
  if ( eps == 0.0L ) return true; // check for both to be zero
  eps *= EPS;
  return ( fabsl( d1 - d2 ) < eps );
} // equal

bool isnotequal( const long double d1, const long double d2 )
// eps hold the minimum distance between the values 
// that will be considered as the numbers are equal
// considering the magnitude of the numbers
{
  const long double EPS = LDBL_EPSILON;
  long double eps1 = fabsl( d1 );
  long double eps2 = fabsl( d2 );
  long double eps = ( eps1 > eps2 ) ? eps1 : eps2;
  if ( eps == 0.0L ) return false; // check for both to be zero
  eps *= EPS;
  return ( fabsl( d1 - d2 ) >= eps );
} // notequal

size_t iamax( const size_t n, const double* const a, const size_t inca )
// return index of maximum absolute value of vector
{
  if ( n < 1 ) return 0;
  if ( n == 1 ) return 1;
  size_t ia = 0;
  double work = fabs( a[ia] );
  size_t index = 1;
  ia += inca;
  for ( size_t i = 2; i <= n; i++ ) {
    if ( fabs( a[ia] ) > work ) {
      work = fabs( a[ia] );
      index = i;
    }
    ia += inca;
  }
  return index;
} // iamax

void swap( const size_t n, double* const a, const size_t inca, double* const b, const size_t incb )
// interchanges two vectors
{
  double s;
  if ( inca != 1 || incb != 1 ) {
    size_t ia = 0;
    size_t ib = 0;
    for ( size_t i = n; i--; ) {
      s = a[ia];
      a[ia] = b[ib];
      b[ib] = s;
      ia += inca;
      ib += incb;
    }
    return;
  }
  size_t i = n >> 2;
  size_t j = 0;
  size_t k = n & 3;
  while ( i-- ) {
    s = a[j]; a[j] = b[j]; b[j] = s;
    s = a[j + 1]; a[j + 1] = b[j + 1]; b[j + 1] = s;
    s = a[j + 2]; a[j + 2] = b[j + 2]; b[j + 2] = s;
    s = a[j + 3]; a[j + 3] = b[j + 3]; b[j + 3] = s;
    j += 4;
  }
  switch ( k ) {
  case 3: s = a[j]; a[j] = b[j]; b[j] = s; j++;
  case 2: s = a[j]; a[j] = b[j]; b[j] = s; j++;
  case 1: s = a[j]; a[j] = b[j]; b[j] = s; j++;
  }
} // swap

double* getvector( const size_t nr, const double c )
// allocates vector space on the heap
{
  double* ptr = 0;
  if ( nr == 0 ) return ptr;
  ptr = ( double* ) calloc( nr, sizeof( double ) );
  ptr--;
  for ( size_t i = 1; i <= nr; i++ ) ptr[i] = c;
  return ptr;
} // getvector

void freevector( double* a )
// de-allocates vector space from the heap
{
  free( ++a );
} // freevector

double** getmatrix( const size_t nr, const size_t nc, const double c )
// allocates matrix space on the heap
{
  double** ptr = 0;
  if ( nr == 0 || nc == 0 ) return ptr;
  double* block = 0;
  ptr = ( double** ) calloc( nr, sizeof( double* ) );
  block = ( double* ) calloc( nr*nc, sizeof( double ) );
  ptr--;
  block--;
  for ( size_t i = 1, im1 = 0; i <= nr; i++, im1++ ) {
    ptr[i] = &block[im1*nc];
    for ( size_t j = 1; j <= nc; j++ ) ptr[i][j] = c;
  }
  return ptr;
} // getmatrix

void freematrix( double** a )
// de-allocates matrix space from the heap
{
  free( ++a[1] ); 
  free( ++a );
} // freematrix

void scal( const size_t n, const long double c, double* const a, const size_t inca )
// scales vector a by constant c
// if c equals 1.0 there is a fast return
// if c equals 0.0 vector a is set to zero
{
  if ( n == 0 ) return;
  if ( isequal( c, 1.0L ) ) return;
  if ( inca != 1 ) {
    size_t ia = 0;
    for ( size_t i = n; i--; ) {
      a[ia] = c * a[ia];
      ia += inca;
    }
    return;
  }
  size_t i = n >> 2;
  size_t j = 0;
  size_t k = n & 3;
  if ( iszero( c ) ) {
    while ( i-- ) {
      a[j] = 0.0;
      a[j + 1] = 0.0;
      a[j + 2] = 0.0;
      a[j + 3] = 0.0;
      j += 4;
    }
    switch ( k ) {
    case 3: a[j] = 0.0; j++;
    case 2: a[j] = 0.0; j++;
    case 1: a[j] = 0.0; j++;
    }
  }
  else {
    register const long double C = c;
    register long double a0, a1, a2, a3;
    while ( i-- ) {
      a0 = a[j];
      a1 = a[j + 1];
      a2 = a[j + 2];
      a3 = a[j + 3];
      a0 *= C;
      a1 *= C;
      a2 *= C;
      a3 *= C;
      a[j] = a0;
      a[j + 1] = a1;
      a[j + 2] = a2;
      a[j + 3] = a3;
      j += 4;
    }
    switch ( k ) {
    case 3: a[j] *= C; j++;
    case 2: a[j] *= C; j++;
    case 1: a[j] *= C; j++;
    }
  }
} // scal

void zeroall( const size_t n, double* const a )
// set all elements of vector a equal to zero
{
  memset( a, 0, n * sizeof( a[1] ) );
} // zeroall

void gemm( const bool transa, const bool transb, const size_t nrc, const size_t ncc, const size_t nab, const long double alpha, double** const a, double** const b, const long double beta, double** const c )
// C = alpha * A(ta) * B(tb) + beta * C
{
  // input cannot be same as output
  assert( a != c );
  assert( b != c );

  // if alpha equals zero
  if ( iszero( alpha ) ) {
    if ( iszero( beta ) ) zeroall( nrc * ncc, &c[1][1] );
    else if ( isnotequal( beta, 1.0L ) ) scal( nrc * ncc, beta, &c[1][1], 1 );
    return;
  }

  if ( isnotzero( beta ) ) scal( nrc * ncc, beta, &c[1][1], 1 );
  else zeroall( nrc * ncc, &c[1][1] );

  // handle transpose options
  if ( transb == false ) {
    if ( transa == false ) {

      // form: C = alpha*A*B + beta*C
      for ( size_t j = 1; j <= ncc; j++ ) {
        for ( size_t k = 1; k <= nab; k++ ) {
          if ( isnotzero( b[k][j] ) ) {
            register long double temp = alpha * b[k][j];
            for ( size_t i = 1; i <= nrc; i++ ) c[i][j] += temp * a[i][k];
          }
        }
      }
    }
    else {

      // form: C = alpha*A'*B + beta*C for beta == 0.0
      if ( iszero( beta ) ) {
        for ( size_t j = 1; j <= ncc; j++ ) {
          for ( size_t i = 1; i <= nrc; i++ ) {
            register long double work = 0.0L;
            for ( size_t k = 1; k <= nab; k++ ) work += a[k][i] * b[k][j];
            c[i][j] = alpha * work;
          }
        }
      }

      // form: C = alpha*A'*B + beta*C for beta != 0.0
      else {
        for ( size_t j = 1; j <= ncc; j++ ) {
          for ( size_t i = 1; i <= nrc; i++ ) {
            register long double work = 0.0L;
            for ( size_t k = 1; k <= nab; k++ ) work += a[k][i] * b[k][j];
            c[i][j] += alpha * work;
          }
        }
      }
    }
  }
  else {
    if ( transa == false ) {

      // form: C = alpha*A*B' + beta*C
      for ( size_t j = 1; j <= ncc; j++ ) {
        for ( size_t k = 1; k <= nab; k++ ) {
          if ( isnotzero( b[j][k] ) ) {
            register long double work = alpha * b[j][k];
            for ( size_t i = 1; i <= nrc; i++ ) c[i][j] += work * a[i][k];
          }
        }
      }
    }
    else {

      // form: C = alpha*A'*B' + beta*C for beta == 0.0
      if ( iszero( beta ) ) {
        for ( size_t j = 1; j <= ncc; j++ ) {
          for ( size_t i = 1; i <= nrc; i++ ) {
            register long double work = 0.0L;
            for ( size_t k = 1; k <= nab; k++ ) work += a[k][i] * b[j][k];
            c[i][j] = alpha * work;
          }
        }
      }

      // form: C = alpha*A'*B' + beta*C for beta != 0.0
      else {
        for ( size_t j = 1; j <= ncc; j++ ) {
          for ( size_t i = 1; i <= nrc; i++ ) {
            register long double work = 0.0L;
            for ( size_t k = 1; k <= nab; k++ ) work += a[k][i] * b[j][k];
            c[i][j] += alpha * work;
          }
        }
      }
    }
  }
} // gemm

void euclidean( const size_t n, const size_t p, double** a, const size_t m, double** b, double** const r )
// compute euclidean distances r between rows of a and b
{
  for ( size_t j = 1; j <= m; j++ ) {
    for ( size_t i = 1; i <= n; i++ ) {
      long double sum = 0.0L;
      for ( size_t k = 1; k <= p; k++ ) {
        long double diff = a[i][k] - b[j][k];
        if ( isnotzero( diff ) ) sum += diff * diff;
      }
      r[i][j] = sqrt( sum );
    }
  }
} // euclidean

int inverse1x1( const size_t n, double** const a )
{
  if ( n != 1 ) return 1;
  long double det = a[1][1];
  if ( iszero( det ) ) return 1;
  a[1][1] = 1.0L / det;
  return 0;
} // inverse1x1

int inverse2x2( const size_t n, double** const a )
{
  if ( n != 2 ) return 1;
  long double det = a[1][1] * a[2][2] - a[1][2] * a[2][1];
  if ( iszero( det ) ) return 1;
  long double temp = a[1][1];
  a[1][1] = a[2][2] / det;
  a[1][2] = -a[1][2] / det;
  a[2][1] = -a[2][1] / det;
  a[2][2] = temp / det;
  return 0;
} // inverse2x2

int inverse3x3( const size_t n, double** const a )
{
  if ( n != 3 ) return 1;
  long double r[9];
  r[0] = a[3][3] * a[2][2] - a[3][2] * a[2][3];
  r[1] = a[3][2] * a[1][3] - a[3][3] * a[1][2];
  r[2] = a[2][3] * a[1][2] - a[2][2] * a[1][3];
  r[3] = a[3][1] * a[2][3] - a[3][3] * a[2][1];
  r[4] = a[3][3] * a[1][1] - a[3][1] * a[1][3];
  r[5] = a[2][1] * a[1][3] - a[2][3] * a[1][1];
  r[6] = a[3][2] * a[2][1] - a[3][1] * a[2][2];
  r[7] = a[3][1] * a[1][2] - a[3][2] * a[1][1];
  r[8] = a[2][2] * a[1][1] - a[2][1] * a[1][2];
  long double det = a[1][1] * r[0] + a[2][1] * r[1] + a[3][1] * r[2];
  if ( iszero( det ) ) return 1;
  det = 1.0L / det;
  for ( size_t i = 1, k = 0; i <= 3; i++ ) for ( size_t j = 1; j <= 3; j++, k++ ) a[i][j] = det * r[k];
  return 0;
} // inverse3x3

int inverse4x4( const size_t n, double** const a )
{
  if ( n != 4 ) return 1;

  long double t[12];
  long double r[16];

  t[0] = a[3][3] * a[4][4];
  t[1] = a[4][3] * a[3][4];
  t[2] = a[2][3] * a[4][4];
  t[3] = a[4][3] * a[2][4];
  t[4] = a[2][3] * a[3][4];
  t[5] = a[3][3] * a[2][4];
  t[6] = a[1][3] * a[4][4];
  t[7] = a[4][3] * a[1][4];
  t[8] = a[1][3] * a[3][4];
  t[9] = a[3][3] * a[1][4];
  t[10] = a[1][3] * a[2][4];
  t[11] = a[2][3] * a[1][4];

  r[0] = t[0] * a[2][2] + t[3] * a[3][2] + t[4] * a[4][2];
  r[0] -= t[1] * a[2][2] + t[2] * a[3][2] + t[5] * a[4][2];
  r[1] = t[1] * a[1][2] + t[6] * a[3][2] + t[9] * a[4][2];
  r[1] -= t[0] * a[1][2] + t[7] * a[3][2] + t[8] * a[4][2];
  r[2] = t[2] * a[1][2] + t[7] * a[2][2] + t[10] * a[4][2];
  r[2] -= t[3] * a[1][2] + t[6] * a[2][2] + t[11] * a[4][2];
  r[3] = t[5] * a[1][2] + t[8] * a[2][2] + t[11] * a[3][2];
  r[3] -= t[4] * a[1][2] + t[9] * a[2][2] + t[10] * a[3][2];
  r[4] = t[1] * a[2][1] + t[2] * a[3][1] + t[5] * a[4][1];
  r[4] -= t[0] * a[2][1] + t[3] * a[3][1] + t[4] * a[4][1];
  r[5] = t[0] * a[1][1] + t[7] * a[3][1] + t[8] * a[4][1];
  r[5] -= t[1] * a[1][1] + t[6] * a[3][1] + t[9] * a[4][1];
  r[6] = t[3] * a[1][1] + t[6] * a[2][1] + t[11] * a[4][1];
  r[6] -= t[2] * a[1][1] + t[7] * a[2][1] + t[10] * a[4][1];
  r[7] = t[4] * a[1][1] + t[9] * a[2][1] + t[10] * a[3][1];
  r[7] -= t[5] * a[1][1] + t[8] * a[2][1] + t[11] * a[3][1];

  t[0] = a[3][1] * a[4][2];
  t[1] = a[4][1] * a[3][2];
  t[2] = a[2][1] * a[4][2];
  t[3] = a[4][1] * a[2][2];
  t[4] = a[2][1] * a[3][2];
  t[5] = a[3][1] * a[2][2];
  t[6] = a[1][1] * a[4][2];
  t[7] = a[4][1] * a[1][2];
  t[8] = a[1][1] * a[3][2];
  t[9] = a[3][1] * a[1][2];
  t[10] = a[1][1] * a[2][2];
  t[11] = a[2][1] * a[1][2];

  r[8] = t[0] * a[2][4] + t[3] * a[3][4] + t[4] * a[4][4];
  r[8] -= t[1] * a[2][4] + t[2] * a[3][4] + t[5] * a[4][4];
  r[9] = t[1] * a[1][4] + t[6] * a[3][4] + t[9] * a[4][4];
  r[9] -= t[0] * a[1][4] + t[7] * a[3][4] + t[8] * a[4][4];
  r[10] = t[2] * a[1][4] + t[7] * a[2][4] + t[10] * a[4][4];
  r[10] -= t[3] * a[1][4] + t[6] * a[2][4] + t[11] * a[4][4];
  r[11] = t[5] * a[1][4] + t[8] * a[2][4] + t[11] * a[3][4];
  r[11] -= t[4] * a[1][4] + t[9] * a[2][4] + t[10] * a[3][4];
  r[12] = t[2] * a[3][3] + t[5] * a[4][3] + t[1] * a[2][3];
  r[12] -= t[4] * a[4][3] + t[0] * a[2][3] + t[3] * a[3][3];
  r[13] = t[8] * a[4][3] + t[0] * a[1][3] + t[7] * a[3][3];
  r[13] -= t[6] * a[3][3] + t[9] * a[4][3] + t[1] * a[1][3];
  r[14] = t[6] * a[2][3] + t[11] * a[4][3] + t[3] * a[1][3];
  r[14] -= t[10] * a[4][3] + t[2] * a[1][3] + t[7] * a[2][3];
  r[15] = t[10] * a[3][3] + t[4] * a[1][3] + t[9] * a[2][3];
  r[15] -= t[8] * a[2][3] + t[11] * a[3][3] + t[5] * a[1][3];

  long double det = a[1][1] * r[0] + a[2][1] * r[1] + a[3][1] * r[2] + a[4][1] * r[3];

  if ( iszero( det ) ) return 1;

  det = 1.0L / det;
  for ( size_t i = 1, k = 0; i <= 4; i++ ) for ( size_t j = 1; j <= 4; j++, k++ ) a[i][j] = det * r[k];

  return 0;
} // inverse4x4

int inverse( const size_t n, double** a )
// compute inverse from real symmetric matrix a
// inverse is returned in matrix a
{
  // speed up computations for small matrices
  if ( n == 1 ) return( inverse1x1( n, a ) );
  if ( n == 2 ) return( inverse2x2( n, a ) );
  if ( n == 3 ) return( inverse3x3( n, a ) );
  if ( n == 4 ) return( inverse4x4( n, a ) );

  size_t* idx = ( size_t* ) calloc( n, sizeof( size_t ) );
  idx--;
  //for ( size_t i = 1; i <= n; i++ ) idx[i] = 0;
  for ( size_t j = 1; j <= n; j++ ) {
    size_t jp = j - 1 + iamax( n - j + 1, &a[j][j], n );
    idx[j] = jp;
    if ( isnotzero( a[jp][j] ) ) {
      if ( jp != j ) swap( n, &a[j][1], 1, &a[jp][1], 1 );
      if ( j < n ) for ( size_t i = j + 1; i <= n; i++ ) a[i][j] /= a[j][j];
    }
    else return j;
    if ( j < n ) {
      for ( size_t i = j + 1; i <= n; i++ ) {
        for ( size_t ii = j + 1; ii <= n; ii++ ) a[i][ii] -= a[i][j] * a[j][ii];
      }
    }
  }
  for ( size_t j = 1; j <= n; j++ ) {
    if ( iszero( a[j][j] ) ) return j;
  }
  for ( size_t j = 1; j <= n; j++ ) {
    a[j][j] = 1.0 / a[j][j];
    long double ajj = -a[j][j];
    for ( size_t i = 1; i <= j - 1; i++ ) {
      if ( isnotequal( a[i][j], 0.0 ) ) {
        long double temp = a[i][j];
        for ( size_t ii = 1; ii <= i - 1; ii++ ) a[ii][j] += temp * a[ii][i];
        a[i][j] *= a[i][i];
      }
    }
    scal( j - 1, ajj, &a[1][j], n );
  }
  double* v = getvector( n, 0.0 );
  for ( size_t j = n - 1; j >= 1; j-- ) {
    for ( size_t i = j + 1; i <= n; i++ ) {
      v[i] = a[i][j];
      a[i][j] = 0.0;
    }
    for ( size_t jj = j + 1, k = 1; jj <= n; jj++, k++ ) {
      long double work = -1.0 * v[jj];
      for ( size_t ii = 1; ii <= n; ii++ ) a[ii][j] += work * a[ii][jj];
    }
  }
  for ( size_t j = n - 1; j >= 1; j-- ) {
    size_t jp = idx[j];
    if ( jp != j ) swap( n, &a[1][j], n, &a[1][jp], n );
  }
  free( ++idx );
  free( ++v );
  return 0;
} // inverse

long double plogis( const long double x )
{ 
  const long double expx = expl( x );
  return( expx / ( 1.0L + expx ) );
} // plogis

void lrmove( const double* const a, double* const b, const size_t n )
// moves vector a to vector b for n elements
{
  size_t i = n >> 3;
  size_t j = 0;
  size_t k = n & 7;
  while ( i-- ) {
    b[j] = a[j];
    b[j + 1] = a[j + 1];
    b[j + 2] = a[j + 2];
    b[j + 3] = a[j + 3];
    b[j + 4] = a[j + 4];
    b[j + 5] = a[j + 5];
    b[j + 6] = a[j + 6];
    b[j + 7] = a[j + 7];
    j += 8;
  }
  switch ( k ) {
  case 7: b[j] = a[j]; j++;
  case 6: b[j] = a[j]; j++;
  case 5: b[j] = a[j]; j++;
  case 4: b[j] = a[j]; j++;
  case 3: b[j] = a[j]; j++;
  case 2: b[j] = a[j]; j++;
  case 1: b[j] = a[j]; j++;
  }
} // lrmove

long double asum( const size_t n, const double* const a, const size_t inca )
// return the sum of the absolute values of vector a
{
  if ( inca != 1 ) {
    size_t ia = 0;
    long double s = 0.0L;
    for ( size_t i = n; i--;) {
      s += fabs( a[ia] );
      ia += inca;
    }
    return s;
  }
  size_t i = n >> 3;
  size_t j = 0;
  size_t k = n & 7;
  long double s = 0.0L;
  while ( i-- ) {
    s += fabs( a[j] ) + fabs( a[j + 1] ) + fabs( a[j + 2] ) + fabs( a[j + 3] ) + fabs( a[j + 4] ) + fabs( a[j + 5] ) + fabs( a[j + 6] ) + fabs( a[j + 7] );
    j += 8;
  }
  switch ( k ) {
  case 7: s += fabs( a[j] ); j++;
  case 6: s += fabs( a[j] ); j++;
  case 5: s += fabs( a[j] ); j++;
  case 4: s += fabs( a[j] ); j++;
  case 3: s += fabs( a[j] ); j++;
  case 2: s += fabs( a[j] ); j++;
  case 1: s += fabs( a[j] ); j++;
  }
  return s;
} // asum

long double dot( const size_t n, const double* const a, const size_t inca, const double* const b, const size_t incb )
// returns the dot product of two vectors
// with the same vector for a and b, the sum-of-squares is returned
{
  if ( n == 0 ) return 0.0L;
  long double s = 0.0L;
  if ( inca != 1 || incb != 1 ) {
    size_t ia = 0;
    size_t ib = 0;
    for ( size_t i = n; i--; ) {
      s += a[ia] * b[ib];
      ia += inca;
      ib += incb;
    }
    return s;
  }
  size_t i = n >> 3;
  size_t j = 0;
  size_t k = n & 7;
  while ( i-- ) {
    s += a[j] * b[j];
    s += a[j + 1] * b[j + 1];
    s += a[j + 2] * b[j + 2];
    s += a[j + 3] * b[j + 3];
    s += a[j + 4] * b[j + 4];
    s += a[j + 5] * b[j + 5];
    s += a[j + 6] * b[j + 6];
    s += a[j + 7] * b[j + 7];
    j += 8;
  }
  switch ( k ) {
  case 7: s += a[j] * b[j]; j++;
  case 6: s += a[j] * b[j]; j++;
  case 5: s += a[j] * b[j]; j++;
  case 4: s += a[j] * b[j]; j++;
  case 3: s += a[j] * b[j]; j++;
  case 2: s += a[j] * b[j]; j++;
  case 1: s += a[j] * b[j]; j++;
  }
  return s;
} // dot

long double sign( const long double r, const long double s )
// returns r with sign of s
{
  if ( s == 0.0L ) return r;
  else return ( s < 0.0L ) ? -r : r;
} // sign

void axpy( const size_t n, const long double c, double* a, const size_t inca, double* b, const size_t incb )
// constant c times vector a plus vector b is returned in vector b: b = b+ca
{
  if ( iszero( c ) ) return;
  if ( inca != 1 || incb != 1 ) {
    size_t ia = 0;
    size_t ib = 0;
    for ( size_t i = n; i--; ) {
      b[ib] += c * a[ia];
      ia += inca;
      ib += incb;
    }
    return;
  }
  size_t i = n >> 3;
  size_t j = 0;
  size_t k = n & 7;
  double* __restrict ra = a;
  double* __restrict rb = b;
  while ( i-- ) {

    const double ca1 = ra[j];
    const double ca2 = ra[j + 1];
    const double ca3 = ra[j + 2];
    const double ca4 = ra[j + 3];
    const double ca5 = ra[j + 4];
    const double ca6 = ra[j + 5];
    const double ca7 = ra[j + 6];
    const double ca8 = ra[j + 7];

    const double cb1 = rb[j];
    const double cb2 = rb[j + 1];
    const double cb3 = rb[j + 2];
    const double cb4 = rb[j + 3];
    const double cb5 = rb[j + 4];
    const double cb6 = rb[j + 5];
    const double cb7 = rb[j + 6];
    const double cb8 = rb[j + 7];

    const double nb1 = cb1 + ( ca1 * c );
    const double nb2 = cb2 + ( ca2 * c );
    const double nb3 = cb3 + ( ca3 * c );
    const double nb4 = cb4 + ( ca4 * c );
    const double nb5 = cb5 + ( ca5 * c );
    const double nb6 = cb6 + ( ca6 * c );
    const double nb7 = cb7 + ( ca7 * c );
    const double nb8 = cb8 + ( ca8 * c );

    rb[j] = nb1;
    rb[j + 1] = nb2;
    rb[j + 2] = nb3;
    rb[j + 3] = nb4;
    rb[j + 4] = nb5;
    rb[j + 5] = nb6;
    rb[j + 6] = nb7;
    rb[j + 7] = nb8;

    j += 8;
  }
  switch ( k ) {
  case 7: b[j] += c * a[j]; j++;
  case 6: b[j] += c * a[j]; j++;
  case 5: b[j] += c * a[j]; j++;
  case 4: b[j] += c * a[j]; j++;
  case 3: b[j] += c * a[j]; j++;
  case 2: b[j] += c * a[j]; j++;
  case 1: b[j] += c * a[j]; j++;
  }
} // axpy

long double pythag( long double x, long double y )
// return sqrt( x * x + y * y ) without problems
{
  long double xabs = fabsl( x );
  long double yabs = fabsl( y );
  long double w = fmaxl( xabs, yabs );
  long double z = fminl( xabs, yabs );
  if ( iszero( z ) ) return( w );
  else {
    long double work = z / w;
    return( w * sqrtl( 1.0L + work * work ) );
  }
} // pythag

int evdcmp( const size_t n, double** eigvec, double* diagnl )
// eigenvalue decomposition eigvec = vwv'
// matrix eigvec is replaced by the eigenvector matrix eigvec
// the eigenvalues are returned in vector diagnl
{
  double* avct = getvector( n * ( n + 1 ) / 2, 0.0 );
  double* subdgl = getvector( n, 0.0 );
  const size_t MAXITER = 40;

  // transfer subdiagonal matrix eigvec to vector avct
  for ( size_t k = 1, i = 1; i <= n; i++ )
    for ( size_t j = 1; j <= i; j++, k++ )
      avct[k] = eigvec[i][j];

  // reduce real symmetric matrix a to a tridiagonal symmetric matrix by householder transformation
  for ( size_t i = n; i >= 1; i-- ) {
    size_t im1 = i - 1;
    size_t ptraij = i * im1 / 2;
    if ( im1 < 1 ) {
      subdgl[1] = 0.0;
      diagnl[1] = avct[1];
      avct[1] = 0.0;
    }
    else {
      lrmove( &avct[ptraij + 1], &diagnl[1], im1 );
      long double scale = asum( im1, &diagnl[1], 1 );
      ptraij += im1;
      if ( iszero( scale ) ) {
        subdgl[i] = 0.0;
        diagnl[i] = avct[ptraij + 1];
        avct[ptraij + 1] = 0.0;
      }
      else {
        scal( im1, 1.0L / scale, &diagnl[1], 1 );
        long double hnum = dot( im1, &diagnl[1], 1, &diagnl[1], 1 );
        long double f = diagnl[im1];
        long double g = -sign( sqrtl( hnum ), f );
        subdgl[i] = scale * g;
        hnum -= f * g;
        diagnl[im1] = f - g;
        avct[ptraij] = scale * diagnl[im1];
        size_t ptri1 = 1;
        for ( size_t k = 1; k <= im1; k++ ) {
          long double dotprd = 0.0L;
          size_t ptrij = ptri1;
          for ( size_t j = 1; j <= k - 1; j++ ) {
            dotprd += avct[ptrij] * diagnl[j];
            ptrij++;
          }
          for ( size_t j = k; j <= im1; j++ ) {
            dotprd += avct[ptrij] * diagnl[j];
            ptrij += j;
          }
          subdgl[k] = dotprd;
          ptri1 += k;
        }
        long double recph = 1.0L / hnum;
        scal( im1, recph, &subdgl[1], 1 );
        long double knum = 0.5L * recph * dot( im1, &subdgl[1], 1, &diagnl[1], 1 );
        size_t ptrarj = 1;
        for ( size_t j = 1; j <= im1; j++ ) {
          f = diagnl[j];
          g = subdgl[j] - knum * f;
          subdgl[j] = g;
          axpy( j, -f, &subdgl[1], 1, &avct[ptrarj], 1 );
          axpy( j, -g, &diagnl[1], 1, &avct[ptrarj], 1 );
          ptrarj += j;
        }
        diagnl[i] = avct[ptraij + 1];
        avct[ptraij + 1] = scale * sqrtl( hnum );
      }
    }
  }

  // set identity matrix eigvec
  for ( size_t i = 1; i <= n; i++ ) {
    for ( size_t j = 1; j <= n; j++ ) eigvec[i][j] = 0.0;
    eigvec[i][i] = 1.0;
  }

  // determine eigenvalues and "eigenvectors" using the ql algorithm with implicit shift
  int retcod = 0;
  if ( n >= 2 ) {
    lrmove( &subdgl[2], &subdgl[1], n - 1 );
    subdgl[n] = 0.0;
  }
  size_t eigenv = 0;
  while ( retcod == 0 && eigenv < n ) {
    eigenv++;
    size_t iter = 0;
    while ( retcod == 0 && iter < MAXITER ) {
      size_t m = 0;
      for ( m = eigenv; m <= n - 1; m++ ) {
        long double work = fabs( diagnl[m] ) + fabs( diagnl[m + 1] );
        if ( fabs( subdgl[m] ) <= ( LDBL_EPSILON * work ) ) break;
      }
      if ( m != eigenv ) {
        iter++;
        long double g = ( diagnl[eigenv + 1] - diagnl[eigenv] ) / ( 2.0L * subdgl[eigenv] );
        long double r = pythag( 1.0L, g );
        g = diagnl[m] - diagnl[eigenv] + subdgl[eigenv] / ( g + sign( r, g ) );
        long double sine = 1.0L;
        long double cosine = 1.0L;
        long double p = 0.0L;
        for ( size_t i = ( m - 1 ); i >= eigenv; i-- ) {
          long double f = sine * subdgl[i];
          long double b = cosine * subdgl[i];
          if ( fabsl( f ) >= fabsl( g ) ) {
            cosine = g / f;
            r = pythag( 1.0L, cosine );
            subdgl[i + 1] = f * r;
            sine = 1.0L / r;
            cosine *= sine;
          }
          else {
            sine = f / g;
            r = pythag( 1.0L, sine );
            subdgl[i + 1] = g * r;
            cosine = 1.0L / r;
            sine *= cosine;
          }
          g = diagnl[i + 1] - p;
          r = ( diagnl[i] - g ) * sine + 2.0L * cosine * b;
          p = sine * r;
          diagnl[i + 1] = g + p;
          g = cosine * r - b;
          for ( size_t k = 1; k <= n; k++ ) {
            b = eigvec[k][i];
            f = eigvec[k][i + 1];
            eigvec[k][i + 1] = sine * b + cosine * f;
            eigvec[k][i] = cosine * b - sine * f;
          }
        }
        diagnl[eigenv] = diagnl[eigenv] - p;
        subdgl[eigenv] = g;
        subdgl[m] = 0.0;
      }
      else retcod = -1;
    }
    if ( retcod == 0 ) retcod = ( int ) ( eigenv );
    else retcod = 0;
  }

  // set eigenvalues in descending order
  // also order corresponding "eigenvectors"
  for ( size_t i = 1; i <= n - 1; i++ ) {
    size_t k = i;
    long double work = diagnl[i];
    for ( size_t j = i + 1; j <= n; j++ ) {
      if ( diagnl[j] > work ) work = diagnl[k = j];
    }
    if ( k != i ) {
      diagnl[k] = diagnl[i];
      diagnl[i] = work;
      for ( size_t j = 1; j <= n; j++ ) {
        long double d = eigvec[j][i];
        eigvec[j][i] = eigvec[j][k];
        eigvec[j][k] = d;
      }
    }
  }

  // transform back the eigenvectors
  size_t ptii = 1;
  for ( size_t i = 2; i <= n; i++ ) {
    size_t pti1 = ptii + 1;
    ptii += i;
    long double aii = avct[ptii];
    if ( isnotzero( aii ) ) {
      size_t im1 = i - 1;
      for ( size_t j = 1; j <= n; j++ ) {
        long double dotprd = dot( im1, &avct[pti1], 1, &eigvec[1][j], n ); // index: j-1
        dotprd = -( dotprd / aii ) / aii;
        axpy( im1, dotprd, &avct[pti1], 1, &eigvec[1][j], n );             // index: j-1
      }
    }
  }

  // free memory and return
  free( ++avct );
  free( ++subdgl );

  return retcod;
} // evdcmp

void set( const size_t n, const double b, double* const a, const size_t inca )
// set elements of vector a equal to scalar b
{
  if ( inca != 1 ) {
    size_t ia = 0;
    for ( size_t i = n; i--; ) {
      a[ia] = b;
      ia += inca;
    }
    return;
  }
  size_t i = n >> 3;
  size_t j = 0;
  size_t k = n & 7;
  while ( i-- ) {
    a[j] = b;
    a[j + 1] = b;
    a[j + 2] = b;
    a[j + 3] = b;
    a[j + 4] = b;
    a[j + 5] = b;
    a[j + 6] = b;
    a[j + 7] = b;
    j += 8;
  }
  switch ( k ) {
    case 7: a[j] = b; j++;
    case 6: a[j] = b; j++;
    case 5: a[j] = b; j++;
    case 4: a[j] = b; j++;
    case 3: a[j] = b; j++;
    case 2: a[j] = b; j++;
    case 1: a[j] = b;
  }
} // set

void copy( const size_t n, const double* const a, const size_t inca, double* const b, const size_t incb )
// copies vector a to vector b
{
  if ( n == 0 ) return;
  if ( inca != 1 || incb != 1 ) {
    size_t ia = 0;
    size_t ib = 0;
    for ( size_t i = n; i--; ) {
      b[ib] = a[ia];
      ia += inca;
      ib += incb;
    }
    return;
  }
  size_t i = n >> 3;
  size_t j = 0;
  size_t k = n & 7;
  while ( i-- ) {
    b[j] = a[j];
    b[j + 1] = a[j + 1];
    b[j + 2] = a[j + 2];
    b[j + 3] = a[j + 3];
    b[j + 4] = a[j + 4];
    b[j + 5] = a[j + 5];
    b[j + 6] = a[j + 6];
    b[j + 7] = a[j + 7];
    j += 8;
  }
  switch ( k ) {
    case 7: b[j] = a[j]; j++;
    case 6: b[j] = a[j]; j++;
    case 5: b[j] = a[j]; j++;
    case 4: b[j] = a[j]; j++;
    case 3: b[j] = a[j]; j++;
    case 2: b[j] = a[j]; j++;
    case 1: b[j] = a[j]; j++;
  }
} // copy

void rotation( const size_t n, const size_t p, double** z, double** r, double* ev )
// return principal axes rotation matrix R
// return identity matrix on error
{
  gemm( true, false, p, p, n, 1.0, z, z, 0.0, r );
  if ( evdcmp( p, r, ev ) != 0 ) {
    set( p * p, 0.0, &r[1][1], 1 );
    for ( size_t k = 1; k <= p; k++ ) r[k][k] = 1.0;
    return;
  }
  for ( size_t k = 1; k <= p; k++ ) {
    long double work = 0.0L;
    for ( size_t i = 1; i <= p; i++ ) work += z[1][i] * r[i][k];
    if ( work < 0.0L ) for ( size_t i = 1; i <= p; i++ ) r[i][k] *= -1.0;
  }
} // rotation

void weightedrotation( const size_t n, const size_t p, double** z, double* w, double** r, double* ev )
// return principal axes rotation matrix R
// return identity matrix on error
{
  for ( size_t i = 1; i <= p; i++ ) {
    for ( size_t j = 1; j <= p; j++ ) {
      long double work = 0.0L;
      for ( size_t k = 1; k <= n; k++ ) work += z[k][i] * w[k] * z[k][j];
      r[i][j] = work;
    }
  }
  if ( evdcmp( p, r, ev ) != 0 ) {
    set( p * p, 0.0, &r[1][1], 1 );
    for ( size_t k = 1; k <= p; k++ ) r[k][k] = 1.0;
    return;
  }
  for ( size_t k = 1; k <= p; k++ ) {
    long double work = 0.0L;
    for ( size_t i = 1; i <= p; i++ ) work += z[1][i] * r[i][k];
    if ( work < 0.0L ) for ( size_t i = 1; i <= p; i++ ) r[i][k] *= -1.0;
  }
} // rotation

void rotateplus( const size_t n, const size_t p, double** z, const size_t n1, double** z1 )
// rotate Z to principal axes, plus
{
  double* ev = getvector( n, 0.0 );
  double** r = getmatrix( p, p, 0.0 );
  rotation( n, p, z, r, ev );
  freevector( ev );
  double** tz = getmatrix( n, p, 0.0 );
  gemm( false, false, n, p, p, 1.0, z, r, 0.0, tz );
  copy( n * p, &tz[1][1], 1, &z[1][1], 1 );
  freematrix( tz );
  double** tz1 = getmatrix( n1, p, 0.0 );
  gemm( false, false, n1, p, p, 1.0, z1, r, 0.0, tz1 );
  copy( n1 * p, &tz1[1][1], 1, &z1[1][1], 1 );
  freematrix( tz1 );
  freematrix( r );
} // rotateplus

void weightedrotateplus( const size_t n, const size_t p, double** z, double* w, const size_t n1, double** z1 )
// rotate Z to principal axes, plus
{
  double* ev = getvector( n, 0.0 );
  double** r = getmatrix( p, p, 0.0 );
  weightedrotation( n, p, z, w, r, ev );
  freevector( ev );
  double** tz = getmatrix( n, p, 0.0 );
  gemm( false, false, n, p, p, 1.0, z, r, 0.0, tz );
  copy( n * p, &tz[1][1], 1, &z[1][1], 1 );
  freematrix( tz );
  double** tz1 = getmatrix( n1, p, 0.0 );
  gemm( false, false, n1, p, p, 1.0, z1, r, 0.0, tz1 );
  copy( n1 * p, &tz1[1][1], 1, &z1[1][1], 1 );
  freematrix( tz1 );
  freematrix( r );
} // rotateplus

void rotateplusplus( const size_t n, const size_t p, double** z, const size_t n1, double** z1, const size_t n2, double** z2 )
// rotate Z to principal axes, plus, plus
{
  double* ev = getvector( n, 0.0 );
  double** r = getmatrix( p, p, 0.0 );
  rotation( n, p, z, r, ev );
  freevector( ev );
  double** tz = getmatrix( n, p, 0.0 );
  gemm( false, false, n, p, p, 1.0, z, r, 0.0, tz );
  copy( n * p, &tz[1][1], 1, &z[1][1], 1 );
  freematrix( tz );
  double** tz1 = getmatrix( n1, p, 0.0 );
  gemm( false, false, n1, p, p, 1.0, z1, r, 0.0, tz1 );
  copy( n1 * p, &tz1[1][1], 1, &z1[1][1], 1 );
  freematrix( tz1 );
  double** tz2 = getmatrix( n2, p, 0.0 );
  gemm( false, false, n2, p, p, 1.0, z2, r, 0.0, tz2 );
  copy( n2 * p, &tz2[1][1], 1, &z2[1][1], 1 );
  freematrix( tz2 );
  freematrix( r );
} // rotateplusplus

void rotateplusplusplus( const size_t n, const size_t p, double** z, const size_t n1, double** z1, const size_t n2, double** z2, const size_t n3, double** z3 )
// rotate Z to principal axes, plus, plus, plus
{
  double* ev = getvector( n, 0.0 );
  double** r = getmatrix( p, p, 0.0 );
  rotation( n, p, z, r, ev );
  freevector( ev );
  double** tz = getmatrix( n, p, 0.0 );
  gemm( false, false, n, p, p, 1.0, z, r, 0.0, tz );
  copy( n * p, &tz[1][1], 1, &z[1][1], 1 );
  freematrix( tz );
  double** tz1 = getmatrix( n1, p, 0.0 );
  gemm( false, false, n1, p, p, 1.0, z1, r, 0.0, tz1 );
  copy( n1 * p, &tz1[1][1], 1, &z1[1][1], 1 );
  freematrix( tz1 );
  double** tz2 = getmatrix( n2, p, 0.0 );
  gemm( false, false, n2, p, p, 1.0, z2, r, 0.0, tz2 );
  copy( n2 * p, &tz2[1][1], 1, &z2[1][1], 1 );
  freematrix( tz2 );
  double** tz3 = getmatrix( n3, p, 0.0 );
  gemm( false, false, n3, p, p, 1.0, z3, r, 0.0, tz3 );
  copy( n3 * p, &tz3[1][1], 1, &z3[1][1], 1 );
  freematrix( tz3 );
  freematrix( r );
} // rotateplusplus

