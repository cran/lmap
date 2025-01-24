
#ifndef LMAP_H
#define LMAP_H

#define R

#include "flib.h"
#include "R.h"

extern double mulnomrowresmduneg( const size_t n, const size_t nc, double** g, const size_t m, const size_t px, double** x, double** b, double** v, double** theta, const size_t MAXINNER, const double FCRIT, const size_t MAXITER, const double DCRIT, size_t* lastiter, double* lastdif );
extern void Cmulnomrowresmduneg( int* rn, int* rnc, double* rg, int* rm, int* rpx, double* rx, double* rb, double* rv, double* rtheta, int* rmaxinner, double* rfcrit, int* rmaxiter, double* rdcrit, double* rdeviance );

extern double mulnomcolresmduneg( const size_t n, const size_t nc, double** g, const size_t m, double** u, const size_t pz, double** z, double**c, double** theta, const size_t MAXINNER, const double FCRIT, const size_t MAXITER, const double DCRIT, size_t* lastiter, double* lastdif );
extern void Cmulnomcolresmduneg( int* rn, int* rnc, double* rg, int* rm, double* ru, int* rpz, double* rz, double* rc, double* rtheta, int* rmaxinner, double* rfcrit, int* rmaxiter, double* rdcrit, double* rdeviance );

extern double mulnomresmduneg( const size_t n, const size_t nc, double** g, const size_t m, const size_t px, double** x, double** b, const size_t pz, double** z, double**c, double** theta, const size_t MAXINNER, const double FCRIT, const size_t MAXITER, const double DCRIT, size_t* lastiter, double* lastdif );
extern void Cmulnomresmduneg( int* rn, int* rnc, double* rg, int* rm, int* rpx, double* rx, double* rb, int* rpz, double* rz, double* rc, double* rtheta, int* rmaxinner, double* rfcrit, int* rmaxiter, double* rdcrit, double* rdeviance );
extern double mulvarbinmduneg( const size_t n, const size_t r, double** y, const size_t m, double** u, double** v, const bool mains, double* mu, const size_t MAXINNER, const double FCRIT, const size_t MAXITER, const double DCRIT, size_t* lastiter, double* lastdif );
extern void Cmulvarbinmduneg( int* rn, int* rr, double* ry, int* rm, double* ru, double* rv, int* rmains, double* rmu, int* rmaxinner, double* rfcrit, int* rmaxiter, double* rdcrit, double* rdeviance );

extern double mulvarbinrowresmduneg( const size_t n, const size_t r, double** y, const size_t p, double** x, const size_t m, double** b, double** v, const bool mains, double* mu, const size_t MAXINNER, const double FCRIT, const size_t MAXITER, const double DCRIT, size_t* lastiter, double* lastdif );
extern void Cmulvarbinrowresmduneg( int* rn, int* rr, double* ry, int* rm, double* ru, int* rp, double* rx, double* rb, int* rmains, double* rmu, int* rmaxinner, double* rfcrit, int* rmaxiter, double* rdcrit, double* rdeviance );

extern double mulvarbincolresmduneg( const size_t n, const size_t r, double** y, const size_t m, double** u, const size_t p, double** x, double** b, const bool mains, double* mu, const size_t MAXINNER, const double FCRIT, const size_t MAXITER, const double DCRIT, size_t* lastiter, double* lastdif );
extern void Cmulvarbincolresmduneg( int* rn, int* rr, double* ry, int* rm, double* ru, int* rp, double* rx, double* rb, int* rmains, double* rmu, int* rmaxiter, double* rdcrit, int* rmaxinner, double* rfcrit, double* rdeviance );

extern double mulvarbinresmduneg( const size_t n, const size_t r, double** y, const size_t m, const size_t pu, double** xu, double** bu, const size_t pv, double** xv, double** bv, const bool mains, double* mu, const size_t MAXINNER, const double FCRIT, const size_t MAXITER, const double DCRIT, size_t* lastiter, double* lastdif );
extern void Cmulvarbinresmduneg( int* rn, int* rr, double* ry, int* rm, int* rpu, double* rxu, double* rbu, int* rpv, double* rxv, double* rbv, int* rmains, double* rmu, int* rmaxiter, double* rdcrit, int* rmaxinner, double* rfcrit, double* rdeviance );

extern double mulvarbinwgtmduneg( const size_t n, const size_t r, double** y, double** w, const size_t m, double** u, double** v, const bool mains, double* mu, const size_t MAXINNER, const double FCRIT, const size_t MAXITER, const double DCRIT, size_t* lastiter, double* lastdif );
extern void Cmulvarbinwgtmduneg( int* rn, int* rr, double* ry, double* rw, int* rm, double* ru, double* rv, int* rmains, double* rmu, int* rmaxiter, double* rdcrit, int* rmaxinner, double* rfcrit, double* rdeviance );

#endif
