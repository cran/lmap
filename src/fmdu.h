//
// Copyright (c) 2020 Frank M.T.A. Busing (e-mail: busing at fsw dot leidenuniv dot nl)
// FreeBSD or 2-Clause BSD or BSD-2 License applies, see Http://www.freebsd.org/copyright/freebsd-license.html
// This is a permissive non-copyleft free software license that is compatible with the GNU GPL.
//

#ifndef FMDU_H
#define FMDU_H

#define R

#include "flib.h"


// mdu.c
extern double mdu( const size_t n, const size_t m, double** delta, const size_t p, double** x, int** fx, double** y, int** fy, double** d, const size_t MAXITER, const double FCRIT, size_t* lastiter, double* lastdif, const bool echo );
extern double wgtmdu( const size_t n, const size_t m, double** delta, double** w, const size_t p, double** x, int** fx, double** y, int** fy, double** d, const size_t MAXITER, const double FCRIT, size_t* lastiter, double* lastdif, const bool echo );

// mduneg.c
extern double mduneg( const size_t n, const size_t m, double** delta, const size_t p, double** x, int** fx, double** y, int** fy, double** d, const size_t MAXITER, const double FCRIT, size_t* lastiter, double* lastdif, const bool echo );
extern double wgtmduneg( const size_t n, const size_t m, double** delta, double** w, const size_t p, double** x, int** fx, double** y, int** fy, double** d, const size_t MAXITER, const double FCRIT, size_t* lastiter, double* lastdif, const bool echo );

extern double rowresmdu( const size_t n, const size_t m, double** delta, const size_t p, const size_t h, double** q, double** b, double** y, int** fy, double** d, const size_t MAXITER, const double FCRIT, size_t* lastiter, double* lastdif, const bool echo );
//extern double penrowresmdu( const size_t n, const size_t m, double** delta, const size_t p, const size_t h, double** q, double** b, double** y, int** fy, double** d, const double rlambda, const double llambda, const double glambda, const size_t MAXITER, const double FCRIT, const double ZCRIT, size_t* lastiter, double* lastdif, const bool echo );
extern double rowreswgtmdu( const size_t n, const size_t m, double** delta, double** w, const size_t p, const size_t h, double** q, double** b, double** y, int** fy, double** d, const size_t MAXITER, const double FCRIT, size_t* lastiter, double* lastdif, const bool echo );

extern double rowresmduneg( const size_t n, const size_t m, double** delta, const size_t p, const size_t h, double** q, double** b, double** y, int** fy, double** d, const size_t MAXITER, const double FCRIT, size_t* lastiter, double* lastdif, const bool echo );
extern double rowreswgtmduneg( const size_t n, const size_t m, double** delta, double** w, const size_t p, const size_t h, double** q, double** b, double** y, int** fy, double** d, const size_t MAXITER, const double FCRIT, size_t* lastiter, double* lastdif, const bool echo );

extern double colresmdu( const size_t n, const size_t m, double** delta, const size_t p, double** x, int** fx, const size_t h, double** q, double** b, double** d, const size_t MAXITER, const double FCRIT, size_t* lastiter, double* lastdif, const bool echo );
//extern double pencolresmdu( const size_t n, const size_t m, double** delta, const size_t p, double** x, int** fx, const size_t h, double** q, double** b, double** d, const double rlambda, const double llambda, const double glambda, const size_t MAXITER, const double FCRIT, const double ZCRIT, size_t* lastiter, double* lastdif, const bool echo );
extern double colreswgtmdu( const size_t n, const size_t m, double** delta, double** w, const size_t p, double** x, int** fx, const size_t h, double** q, double** b, double** d, const size_t MAXITER, const double FCRIT, size_t* lastiter, double* lastdif, const bool echo );

extern double colresmduneg( const size_t n, const size_t m, double** delta, const size_t p, double** x, int** fx, const size_t h, double** q, double** b, double** d, const size_t MAXITER, const double FCRIT, size_t* lastiter, double* lastdif, const bool echo );
extern double colreswgtmduneg( const size_t n, const size_t m, double** delta, double** w, const size_t p, double** x, int** fx, const size_t h, double** q, double** b, double** d, const size_t MAXITER, const double FCRIT, size_t* lastiter, double* lastdif, const bool echo );

extern double resmdu( const size_t n, const size_t m, double** delta, const size_t p, const size_t hx, double** qx, double** bx, const size_t hy, double** qy, double** by, double** d, const size_t MAXITER, const double FCRIT, size_t* lastiter, double* lastdif, const bool echo );
extern double reswgtmdu( const size_t n, const size_t m, double** delta, double** w, const size_t p, const size_t hx, double** qx, double** bx, const size_t hy, double** qy, double** by, double** d, const size_t MAXITER, const double FCRIT, size_t* lastiter, double* lastdif, const bool echo );

extern double resmduneg( const size_t n, const size_t m, double** delta, const size_t p, const size_t hx, double** qx, double** bx, const size_t hy, double** qy, double** by, double** d, const size_t MAXITER, const double FCRIT, size_t* lastiter, double* lastdif, const bool echo );
extern double reswgtmduneg( const size_t n, const size_t m, double** delta, double** w, const size_t p, const size_t hx, double** qx, double** bx, const size_t hy, double** qy, double** by, double** d, const size_t MAXITER, const double FCRIT, size_t* lastiter, double* lastdif, const bool echo );

extern double external( const size_t n, const size_t m, double** delta, double** w, const size_t p, double** const fixed, double** const z, double** d, const size_t MAXITER, const double FCRIT, size_t* lastiter, double* lastdif, const bool echo );

extern void CRultrafastmdu( int* rn, int* rm, double* rdata, int* rp, double* rx, double* ry, int* rnsteps, double* rminrate, int* rseed );
extern void CRultrafastmdufxd( int* rn, int* rm, double* rdata, int* rp, double* rx, int* rfx, double* ry, int* rfy, int* rnsteps, double* rminrate, int* rseed );
extern void CRultrafastwgtmdu( int* rn, int* rm, double* rdata, double* rw, int* rp, double* rx, double* ry, int* rnsteps, double* rminrate, int* rseed );
extern void CRultrafastwgtmdufxd( int* rn, int* rm, double* rdata, double* rw, int* rp, double* rx, int* rfx, double* ry, int* rfy, int* rnsteps, double* rminrate, int* rseed );
extern void CRultrafastrowresmdu( int* rn, int* rm, double* rdata, int* rp, int* rh, double* rq, double* rb, double* ry, int* rnsteps, double* rminrate, int* rseed );

#endif
