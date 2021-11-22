#ifndef CIRCLE_H
#define CIRCLE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "debugmacros.h"
#include "double2.h"

#ifndef CIRCLE_NUM_CLOSED_CURVES
#define CIRCLE_NUM_CLOSED_CURVES 1
#endif

#ifndef CIRCLE_RADIUS
#define CIRCLE_RADIUS 1
#endif

#ifndef CIRCLE_SPLINE_DEGREE
#define CIRCLE_SPLINE_DEGREE 2
#endif

#ifndef CIRCLE_NUM_CONTROL_POINTS
#define CIRCLE_NUM_CONTROL_POINTS 50
#endif

#if (CIRCLE_SPLINE_DEGREE != 2 && CIRCLE_SPLINE_DEGREE != 3 && CIRCLE_SPLINE_DEGREE != 4)
#error Undefined.
#endif

#if (CIRCLE_NUM_CONTROL_POINTS < 0)
#error Out of domain.
#endif

#if defined(CIRCLE) // solve a test problem
#if (CIRCLE_NUM_CLOSED_CURVES != 1)
#error CIRCLE_NUM_CLOSED_CURVES must be 1 when CIRCLE is enabled
#endif
#if (CIRCLE_RADIUS != 1)
#error CIRCLE_RADIUS must be 1 when CIRCLE is enabled
#endif
#endif

#ifdef __cplusplus
extern "C" {
#endif
  void isolap2d_input_circle(Isolap2d *p);
  void isolap2d_input_circle_info();
#ifdef __cplusplus
}
#endif

#endif /* CIRCLE_H */
