#ifndef CSHAPE_H
#define CSHAPE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "debugmacros.h"
#include "double2.h"

#ifndef CSHAPE_SPLINE_DEGREE
#define CSHAPE_SPLINE_DEGREE 2
#endif

#ifndef CSHAPE_LATTICE
#define CSHAPE_LATTICE 1.0
#endif

#ifndef CSHAPE_SCALE
#define CSHAPE_SCALE 0.8
#endif

#ifndef CSHAPE_NUM_FOR_X
#define CSHAPE_NUM_FOR_X 5
#endif

#ifndef CSHAPE_NUM_FOR_Y
#define CSHAPE_NUM_FOR_Y 5
#endif

#ifndef CSHAPE_NUM_PER_UNIT
#define CSHAPE_NUM_PER_UNIT 10
#endif

#ifdef __cplusplus
extern "C" {
#endif
  void isolap2d_input_cshape(Isolap2d *p);
  void isolap2d_input_cshape_info();
#ifdef __cplusplus
}
#endif

#endif /* CSHAPE_H */
