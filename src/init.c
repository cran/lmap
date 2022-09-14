#include <stdlib.h>
#include <R_ext/Rdynload.h>


extern void Cmulnomrowresmduneg( int* rn, int* rc, double* rg, int* rp, double* rx, int* rm, double* rb, double* rv, double* ru, double* rtheta, int* rmaxiter, double* rdcrit, int* rmaxinner, double* rfcrit, double* rdeviance );
extern void Cmulvarbinmduneg( int* rn, int* rr, double* ry, int* rm, double* ru, double* rv, int* rmains, double* rmu, int* rmaxiter, double* rdcrit, int* rmaxinner, double* rfcrit, double* rdeviance );
extern void Cmulvarbinrowresmduneg( int* rn, int* rr, double* ry, int* rp, double* rx, int* rm, double* rb, double* rv, int* rmains, double* rmu, int* rmaxiter, double* rdcrit, int* rmaxinner, double* rfcrit, double* rdeviance );
extern void Cmulvarbincolresmduneg( int* rn, int* rr, double* ry, int* rm, double* ru, int* rp, double* rx, double* rb, int* rmains, double* rmu, int* rmaxiter, double* rdcrit, int* rmaxinner, double* rfcrit, double* rdeviance );
extern void Cmulvarbinresmduneg( int* rn, int* rr, double* ry, int* rm, int* rpu, double* rxu, double* rbu, int* rpv, double* rxv, double* rbv, int* rmains, double* rmu, int* rmaxiter, double* rdcrit, int* rmaxinner, double* rfcrit, double* rdeviance );

extern void Ccolresmduneg( int* rn, int* rm, double* rdelta, int* rp, double* rx, int* rh, double* rq, double* rb, double* rd, int* rmaxiter, double* rfdif, double* rfvalue );
extern void Cmduneg( int* rn, int* rm, double* rdelta, int* rp, double* rx, double* ry, double* rd, int* rmaxiter, double* rfdif, double* rfvalue );
extern void Cresmduneg( int* rn, int* rm, double* rdelta, int* rp, int* rhx, double* rqx, double* rbx, int* rhy, double* rqy, double* rby, double* rd, int* rmaxiter, double* rfdif, double* rfvalue );
extern void Crowresmdu( int* rn, int* rm, double* rdelta, int* rp, int* rh, double* rq, double* rb, double* ry, double* rd, int* rmaxiter, double* rfdif, double* rfvalue );
extern void Crowresmduneg( int* rn, int* rm, double* rdelta, int* rp, int* rh, double* rq, double* rb, double* ry, double* rd, int* rmaxiter, double* rfdif, double* rfvalue );

extern void Cwgtmduneg( int* rn, int* rm, double* rdelta, double* rw, int* rp, double* rx, double* ry, double* rd, int* rmaxiter, double* rfdif, double* rfvalue );
extern void Cmulvarbinwgtmduneg( int* rn, int* rr, double* ry, double* rw, int* rm, double* ru, double* rv, int* rmains, double* rmu, int* rmaxiter, double* rdcrit, int* rmaxinner, double* rfcrit, double* rdeviance );
extern void Cmdu( int* rn, int* rm, double* rdelta, int* rp, double* rx, double* ry, double* rd, int* rmaxiter, double* rfdif, double* rfvalue );

extern void Ccolreswgtmduneg( int* rn, int* rm, double* rdelta, double* rw, int* rp, double* rx, int* rh, double* rq, double* rb, double* rd, int* rmaxiter, double* rfdif, double* rfvalue );



static const R_CMethodDef CEntries[] = {
  {"Cmulnomrowresmduneg",      ( DL_FUNC ) &Cmulnomrowresmduneg,         15},
  {"Cmulvarbinmduneg",      ( DL_FUNC ) &Cmulvarbinmduneg,         13},
  {"Cmulvarbinrowresmduneg",      ( DL_FUNC ) &Cmulvarbinrowresmduneg,         15},
  {"Cmulvarbincolresmduneg",      ( DL_FUNC ) &Cmulvarbincolresmduneg,         15},
  {"Cmulvarbinresmduneg",      ( DL_FUNC ) &Cmulvarbinresmduneg,         17},
  {"Ccolresmduneg",      ( DL_FUNC ) &Ccolresmduneg,         12},
  {"Cmduneg",      ( DL_FUNC ) &Cmduneg,         10},
  {"Cresmduneg",      ( DL_FUNC ) &Cresmduneg,         14},
  {"Crowresmdu",      ( DL_FUNC ) &Crowresmdu,         12},
  {"Crowresmduneg",      ( DL_FUNC ) &Crowresmduneg,         12},
  {"Cwgtmduneg",      ( DL_FUNC ) &Cwgtmduneg,         11},
  {"Cmulvarbinwgtmduneg",      ( DL_FUNC ) &Cmulvarbinwgtmduneg,         14},
  {"Cmdu",      ( DL_FUNC ) &Cmdu,         10},
  {"Ccolreswgtmduneg",      ( DL_FUNC ) &Ccolreswgtmduneg,         13},
  {NULL, NULL, 0}
};




void R_init_lmap( DllInfo *dll )
{
  R_registerRoutines( dll, CEntries, NULL, NULL, NULL );
  R_useDynamicSymbols( dll, FALSE );
}
