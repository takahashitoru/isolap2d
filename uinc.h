#ifndef UINC_H
#define UINC_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "double2.h"
#include "mathmacros.h"
#include "debugmacros.h"
#include "isolap2d_defs.h"

#define NULL_UINC ( 0 )

#ifdef __cplusplus
extern "C" {
#endif

  void mkuinc(const int ncol, const int ninc, const double2 *xc, double *bvec);
  void isolap2d_mkuinc(const Isolap2d *p, const int ncol, double *bvec);

#ifdef __cplusplus
}
#endif

#endif /* UINC_H */
