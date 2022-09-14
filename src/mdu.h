
#ifndef MDU_H
#define MDU_H

#include "R.h"

extern long double mdu( const size_t n, const size_t m, double** delta, const size_t p, double** x, double** y, double** d, const size_t MAXITER, const long double FCRIT, size_t* lastiter, long double* lastdif );
extern long double wgtmdu( const size_t n, const size_t m, double** delta, double** w, const size_t p, double** x, double** y, double** d, const size_t MAXITER, const long double FCRIT, size_t* lastiter, long double* lastdif );

extern long double mduneg( const size_t n, const size_t m, double** delta, const size_t p, double** x, double** y, double** d, const size_t MAXITER, const long double FCRIT, size_t* lastiter, long double* lastdif );
extern long double wgtmduneg( const size_t n, const size_t m, double** delta, double** w, const size_t p, double** x, double** y, double** d, const size_t MAXITER, const long double FCRIT, size_t* lastiter, long double* lastdif );

extern long double rowresmdu( const size_t n, const size_t m, double** delta, const size_t p, const size_t h, double** q, double** b, double** y, double** d, const size_t MAXITER, const long double FCRIT, size_t* lastiter, long double* lastdif );
extern long double rowreswgtmdu( const size_t n, const size_t m, double** delta, double** w, const size_t p, const size_t h, double** q, double** b, double** y, double** d, const size_t MAXITER, const long double FCRIT, size_t* lastiter, long double* lastdif );

extern long double rowresmduneg( const size_t n, const size_t m, double** delta, const size_t p, const size_t h, double** q, double** b, double** y, double** d, const size_t MAXITER, const long double FCRIT, size_t* lastiter, long double* lastdif );
extern long double rowreswgtmduneg( const size_t n, const size_t m, double** delta, double** w, const size_t p, const size_t h, double** q, double** b, double** y, double** d, const size_t MAXITER, const long double FCRIT, size_t* lastiter, long double* lastdif );

extern long double colresmdu( const size_t n, const size_t m, double** delta, const size_t p, double** x, const size_t h, double** q, double** b, double** d, const size_t MAXITER, const long double FCRIT, size_t* lastiter, long double* lastdif );
extern long double colreswgtmdu( const size_t n, const size_t m, double** delta, double** w, const size_t p, double** x, const size_t h, double** q, double** b, double** d, const size_t MAXITER, const long double FCRIT, size_t* lastiter, long double* lastdif );

extern long double colresmduneg( const size_t n, const size_t m, double** delta, const size_t p, double** x, const size_t h, double** q, double** b, double** d, const size_t MAXITER, const long double FCRIT, size_t* lastiter, long double* lastdif );
extern long double colreswgtmduneg( const size_t n, const size_t m, double** delta, double** w, const size_t p, double** x, const size_t h, double** q, double** b, double** d, const size_t MAXITER, const long double FCRIT, size_t* lastiter, long double* lastdif );

extern long double resmdu( const size_t n, const size_t m, double** delta, const size_t p, const size_t hx, double** qx, double** bx, const size_t hy, double** qy, double** by, double** d, const size_t MAXITER, const long double FCRIT, size_t* lastiter, long double* lastdif );
extern long double reswgtmdu( const size_t n, const size_t m, double** delta, double** w, const size_t p, const size_t hx, double** qx, double** bx, const size_t hy, double** qy, double** by, double** d, const size_t MAXITER, const long double FCRIT, size_t* lastiter, long double* lastdif );

extern long double resmduneg( const size_t n, const size_t m, double** delta, const size_t p, const size_t hx, double** qx, double** bx, const size_t hy, double** qy, double** by, double** d, const size_t MAXITER, const long double FCRIT, size_t* lastiter, long double* lastdif );
extern long double reswgtmduneg( const size_t n, const size_t m, double** delta, double** w, const size_t p, const size_t hx, double** qx, double** bx, const size_t hy, double** qy, double** by, double** d, const size_t MAXITER, const long double FCRIT, size_t* lastiter, long double* lastdif );

#endif
