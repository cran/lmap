
#ifndef LIBRARY_H
#define LIBRARY_H

#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define iszero( x ) ( ( ( x ) < LDBL_EPSILON ) && ( ( x ) > -LDBL_EPSILON ) )
#define isnotzero( x ) ( ( ( x ) > LDBL_EPSILON ) || ( ( x ) < -LDBL_EPSILON ) )

extern bool isequal( const long double d1, const long double d2 );
extern bool isnotequal( const long double d1, const long double d2 );
extern size_t iamax( const size_t n, const double* const a, const size_t inca );
extern void swap( const size_t n, double* const a, const size_t inca, double* const b, const size_t incb );
extern double* getvector( const size_t nr, const double c );
extern void freevector( double* a );
extern double** getmatrix( const size_t nr, const size_t nc, const double c );
extern void freematrix( double** a );
extern void scal( const size_t n, const long double c, double* const a, const size_t inca );
extern void zeroall( const size_t n, double* const a );
extern void gemm( const bool transa, const bool transb, const size_t nrc, const size_t ncc, const size_t nab, const long double alpha, double** const a, double** const b, const long double beta, double** const c );
extern void euclidean( const size_t n, const size_t p, double** a, const size_t m, double** b, double** const r );
extern int inverse1x1( const size_t n, double** const a );
extern int inverse2x2( const size_t n, double** const a );
extern int inverse3x3( const size_t n, double** const a );
extern int inverse4x4( const size_t n, double** const a );
extern int inverse( const size_t n, double** a );
extern long double plogis( const long double x );
extern void lrmove( const double* const a, double* const b, const size_t n );
extern long double asum( const size_t n, const double* const a, const size_t inca );
extern long double dot( const size_t n, const double* const a, const size_t inca, const double* const b, const size_t incb );
extern long double sign( const long double r, const long double s );
extern void axpy( const size_t n, const long double c, double* a, const size_t inca, double* b, const size_t incb );
extern long double pythag( long double x, long double y );
extern int evdcmp( const size_t n, double** eigvec, double* diagnl );
extern void rotateplus( const size_t n, const size_t p, double** z, const size_t n1, double** z1 );
extern void weightedrotateplus( const size_t n, const size_t p, double** z, double* w, const size_t n1, double** z1 );
extern void rotateplusplus( const size_t n, const size_t p, double** z, const size_t n1, double** z1, const size_t n2, double** z2 );
extern void rotateplusplusplus( const size_t n, const size_t p, double** z, const size_t n1, double** z1, const size_t n2, double** z2, const size_t n3, double** z3 );

#endif
