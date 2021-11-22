#ifndef GAUSSONE_H
#define GAUSSONE_H

#include <stdio.h>
#include <stdlib.h>

#include "debugmacros.h"

#define NO_CENTER (- 1)
#if (NO_CENTER >= 0)
#error NO_CENTER must be less than 0.
#endif

#define GAUSSONE_NMAX ( 24 )

typedef struct {
  int n;
  int c;
  double *x;
  double *w;
} GaussOne;

#ifdef __cplusplus
extern "C" {
#endif
  void GaussOne_load(const int n, GaussOne **g);
  void GaussOne_check(const GaussOne *g);
  void GaussOne_free(GaussOne *g);
#ifdef __cplusplus
}
#endif

#endif /* GAUSSONE_H */
