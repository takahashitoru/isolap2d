#ifndef CLOSEDCURVE_H
#define CLOSEDCURVE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "double2.h"
#include "mathmacros.h"
#include "debugmacros.h"

#if(0) // probably okay, though results will differ by round-off

#define POLY1(x, c0, c1) ( (c0) + (x) * (c1) )
#define POLY2(x, c0, c1, c2) ( (c0) + (x) * ((c1) + (x) * (c2)) )
#define POLY3(x, c0, c1, c2, c3) ( (c0) + (x) * ((c1) + (x) * ((c2) + (x) * (c3))) )
#define POLY4(x, c0, c1, c2, c3, c4) ( (c0) + (x) * ((c1) + (x) * ((c2) + (x) * ((c3) + (x) * ((c4) + (x))))) )

/* B-spline for p=2 (quadratic) */
#define N20(x) ( POLY2(x,   0,   0,   1) / 2 )
#define N21(x) ( POLY2(x, - 3,   6, - 2) / 2 )
#define N22(x) ( POLY2(x,   9, - 6,   1) / 2 )

#define DN20(x) ( POLY1(x,   0,   1) )
#define DN21(x) ( POLY1(x,   3, - 2) )
#define DN22(x) ( POLY1(x, - 3,   1) )

#define DDN20(x) (   1 )
#define DDN21(x) ( - 2 )
#define DDN22(x) (   1 )

#define DDDN20(x) ( 0 )
#define DDDN21(x) ( 0 )
#define DDDN22(x) ( 0 )

/* B-spline for p=3 (cubic) */
#define N30(x) ( POLY3(x,    0,     0,    0,    1) / 6 )
#define N31(x) ( POLY3(x,    4,  - 12,   12,  - 3) / 6 )
#define N32(x) ( POLY3(x, - 44,    60, - 24,    3) / 6 )
#define N33(x) ( POLY3(x,   64,  - 48,   12,  - 1) / 6 )

#define DN30(x) ( POLY2(x,    0,    0,    1) / 2 )
#define DN31(x) ( POLY2(x,  - 4,    8,  - 3) / 2 )
#define DN32(x) ( POLY2(x,   20, - 16,    3) / 2 )
#define DN33(x) ( POLY2(x, - 16,    8,  - 1) / 2 )

#define DDN30(x) ( POLY1(x,   0,   1) )
#define DDN31(x) ( POLY1(x,   4, - 3) )
#define DDN32(x) ( POLY1(x, - 8,   3) )
#define DDN33(x) ( POLY1(x,   4, - 1) )

#define DDDN30(x) (   1 )
#define DDDN31(x) ( - 3 )
#define DDDN32(x) (   3 )
#define DDDN33(x) ( - 1 )

/* B-spline for p=4 (quatic; thanks to Y. Nakai) */
#define N40(x) ( POLY4(x,     0,     0,     0,     0,     1) / 24 )
#define N41(x) ( POLY4(x,   - 5,    20,  - 30,    20,   - 4) / 24 )
#define N42(x) ( POLY4(x,   155, - 300,   210,  - 60,     6) / 24 )
#define N43(x) ( POLY4(x, - 655,   780, - 330,    60,   - 4) / 24 )
#define N44(x) ( POLY4(x,   625, - 500,   150,  - 20,     1) / 24 )

#define DN40(x) ( POLY3(x,    0,     0,     0,     1) / 6 )
#define DN41(x) ( POLY3(x,    5,  - 15,    15,   - 4) / 6 )
#define DN42(x) ( POLY3(x, - 75,   105,  - 45,     6) / 6 )
#define DN43(x) ( POLY3(x,  195, - 165,    45,   - 4) / 6 )
#define DN44(x) ( POLY3(x, - 75,    75,  - 15,     1) / 6 )

#define DDN40(x) ( POLY2(x,    0,    0,   1) / 2 )
#define DDN41(x) ( POLY2(x,  - 5,   10, - 4) / 2 )
#define DDN42(x) ( POLY2(x,   35, - 30,   6) / 2 )
#define DDN43(x) ( POLY2(x, - 55,   30, - 4) / 2 )
#define DDN44(x) ( POLY2(x,   25, - 10,   1) / 2 )

#define DDDN40(x) ( POLY1(x,    0,    1) )
#define DDDN41(x) ( POLY1(x,    5,  - 4) )
#define DDDN42(x) ( POLY1(x, - 15,    6) )
#define DDDN43(x) ( POLY1(x,   15,  - 4) )
#define DDDN44(x) ( POLY1(x,  - 5,    1) )

#endif

typedef struct {
  int p;       // degree of B-spline functions; two or more
  int n;       // number of control points
  double *t;   // knots; t[0:n+p]
  double2 *cp; // control points; cp[0:n-1], where cp[0]=cp[n-p], cp[1]=cp[n-p+1],..., cp[p-1]=cp[n-1]
  double2 x;   // centre or representative point of this closed curve (this is optional)
} Closedcurve;

#ifdef __cplusplus
extern "C" {
#endif

  Closedcurve *closedcurve_allocate(const int ncc);
  void closedcurve_free(Closedcurve *cc);
  double2 *closedcurve_cp_alloc(const int n);
  void closedcurve_cp_free(double2 *cp);

  double *closedcurve_knot_alloc(const int p, const int n);
  void closedcurve_knot_free(double *t);
  void closedcurve_knot_comp(const int p, const int n, double *t);
  void closedcurve_knot_check(const int ncc, const Closedcurve *cc);
  void closedcurve_draw_boundary(const int ncc, const Closedcurve *cc);

  double closedcurve_basis(const int p, const int n, const double *t, const int j, const double tnow);
  double closedcurve_basis_derivative(const int p, const int n, const double *t, const int j, const double tnow);
  double closedcurve_potential(const int p, const int n, const double *t, const double *xvec, const double tnow);
  double2 closedcurve_position(const int p, const int n, const double *t, const double2 *cp, const double tnow);
  double2 closedcurve_position_derivative(const int p, const int n, const double *t, const double2 *cp, const double tnow);
  double2 closedcurve_position_derivative2(const int p, const int n, const double *t, const double2 *cp, const double tnow);
  double2 closedcurve_position_derivative3(const int p, const int n, const double *t, const double2 *cp, const double tnow);

#ifdef __cplusplus
}
#endif


#endif /* CLOSEDCURVE_H */
