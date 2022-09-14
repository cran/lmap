
#ifndef LMDU_H
#define LMDU_H

#include "R.h"

extern void Cmulnomrowresmduneg( int* rn, int* rc, double* rg, int* rp, double* rx, int* rm, double* rb, double* rv, double* ru, double* rtheta, int* rmaxiter, double* rdcrit, int* rmaxinner, double* rfcrit, double* rdeviance );

extern void Cmulvarbinmduneg( int* rn, int* rr, double* ry, int* rm, double* ru, double* rv, int* rmains, double* rmu, int* rmaxiter, double* rdcrit, int* rmaxinner, double* rfcrit, double* rdeviance );
extern void Cmulvarbinrowresmduneg( int* rn, int* rr, double* ry, int* rm, double* ru, int* rp, double* rx, double* rb, int* rmains, double* rmu, int* rmaxiter, double* rdcrit, int* rmaxinner, double* rfcrit, double* rdeviance );
extern void Cmulvarbincolresmduneg( int* rn, int* rr, double* ry, int* rm, double* ru, int* rp, double* rx, double* rb, int* rmains, double* rmu, int* rmaxiter, double* rdcrit, int* rmaxinner, double* rfcrit, double* rdeviance );
extern void Cmulvarbinresmduneg( int* rn, int* rr, double* ry, int* rm, int* rpu, double* rxu, double* rbu, int* rpv, double* rxv, double* rbv, int* rmains, double* rmu, int* rmaxiter, double* rdcrit, int* rmaxinner, double* rfcrit, double* rdeviance );

extern void Cmulvarbinwgtmduneg( int* rn, int* rr, double* ry, double* rw, int* rm, double* ru, double* rv, int* rmains, double* rmu, int* rmaxiter, double* rdcrit, int* rmaxinner, double* rfcrit, double* rdeviance );

#endif
