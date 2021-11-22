#ifndef MSOLVE_H
#define MSOLVE_H

#include <stdio.h>
#include <stdlib.h>

#include "clap2d_defs.h" // defines Clap2d
#include "double2.h" // defines double2
#include "ldedat.h"
#include "direct.h"
#include "bval.h"

#ifdef __cplusplus
extern "C" {
#endif

  void clap2d_msolve(const int n, const void *M, const double *r, double *z);

#ifdef __cplusplus
}
#endif

#endif /* MSOLVE_H */
