#ifndef DXF_H
#define DXF_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "double2.h"
#include "debugmacros.h"

#define DXF_MAX_COLUMNS 256
#define DXF_TRUE 1
#define DXF_FALSE 0

#ifndef DXF_MAX_SPLINES
#define DXF_MAX_SPLINES 10000
#endif

typedef struct {
  int degree; // degree of B-spline functions; two or more
  int flag;   // flag of B-splines
  int *t;     // knots (ignored)
  int m;      // number of control points (no wrapping)
  double *x;  // x-coordinates of control points
  double *y;  // y-coordinates of control points
  double *z;  // z-coordinates of control points
  int cx;     // counter of x-coordinates
  int cy;     // counter of y-coordinates
  int cz;     // counter of z-coordinates
} Bspline;


#ifdef __cplusplus
extern "C" {
#endif
  //120122  void read_bspline_from_dxf(const char *dxffile, int *p, int *n, double2 **cp);
  void dxf_read_bspline(const char *dxffile, const int maxbsp, int *nbsp, Bspline *bsp);
  void dxf_free_bspline(const int nbsp, Bspline *bsp);
#ifdef __cplusplus
}
#endif


#endif /* DXF_H */
