#ifndef HEXA_H
#define HEXA_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "debugmacros.h"
#include "double2.h"

#define SQRT3 (sqrt(3.0))

#ifndef HEXA_SPLINE_DEGREE
#define HEXA_SPLINE_DEGREE 2
#endif

#ifndef HEXA_NUM_CONTROL_POINTS_PER_EDGE
#define HEXA_NUM_CONTROL_POINTS_PER_EDGE 10
#endif

#if (HEXA_NUM_CONTROL_POINTS_PER_EDGE < HEXA_SPLINE_DEGREE)
#error HEXA_NUM_CONTROL_POINTS_PER_EDGE >= HEXA_SPLINE_DEGREE.
#endif

#ifndef HEXA_RADIUS
#define HEXA_RADIUS 1.0
#endif

#ifndef HEXA_SCALE
#define HEXA_SCALE 0.8
#endif

#ifndef HEXA_NUM_FOR_X
#define HEXA_NUM_FOR_X 5
#endif

#ifndef HEXA_NUM_FOR_Y
#define HEXA_NUM_FOR_Y 5
#endif

#ifdef __cplusplus
extern "C" {
#endif
  void isolap2d_input_hexa(Isolap2d *p);
  void isolap2d_input_hexa_info();
#ifdef __cplusplus
}
#endif

#endif /* HEXA_H */
