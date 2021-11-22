#ifndef SSHAPE_H
#define SSHAPE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "debugmacros.h"
#include "double2.h"

#ifndef SSHAPE_SPLINE_DEGREE
#define SSHAPE_SPLINE_DEGREE 2
#endif

#ifndef SSHAPE_LATTICE
#define SSHAPE_LATTICE (1.0)
#endif
#ifndef SSHAPE_SCALE
#define SSHAPE_SCALE (0.95)
#endif
#ifndef SSHAPE_WIDTH
#define SSHAPE_WIDTH (0.05)
#endif
#ifndef SSHAPE_NDIV
#define SSHAPE_NDIV 10
#endif

#ifndef SSHAPE_XRAD
#define SSHAPE_XRAD (SSHAPE_LATTICE * 0.50)
#endif
#ifndef SSHAPE_YRAD
#define SSHAPE_YRAD (SSHAPE_LATTICE * 0.25)
#endif

#ifndef SSHAPE_TOP_STA_ANGLE
#define SSHAPE_TOP_STA_ANGLE (- 0.25 * PI)
#endif
#ifndef SSHAPE_TOP_END_ANGLE
#define SSHAPE_TOP_END_ANGLE (1.5 * PI)
#endif
#ifndef SSHAPE_BOT_STA_ANGLE
#define SSHAPE_BOT_STA_ANGLE (0.5 * PI)
#endif
#ifndef SSHAPE_BOT_END_ANGLE
#define SSHAPE_BOT_END_ANGLE (- 1.25 * PI)
#endif

#ifndef SSHAPE_NUM_FOR_X
#define SSHAPE_NUM_FOR_X 1
#endif
#ifndef SSHAPE_NUM_FOR_Y
#define SSHAPE_NUM_FOR_Y 1
#endif


#ifdef __cplusplus
extern "C" {
#endif
  void isolap2d_input_sshape(Isolap2d *p);
  void isolap2d_input_sshape_info();
#ifdef __cplusplus
}
#endif

#endif /* CSHAPE_H */
