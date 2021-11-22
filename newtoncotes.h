#ifndef NEWTONCOTES_H
#define NEWTONCOTES_H

#include <stdio.h>
#include <stdlib.h>

#include "debugmacros.h"

typedef struct {
  int n; // number of the points, i.e., x_0, x_1, ..., x_{n-1}
  //  double *x; // abcissa
  double h; // equally space of abscissas, i.e., h=1/(n-1)
  double *w; // weight
} Newtoncotes;

#ifdef __cplusplus
extern "C" {
#endif
  void newtoncotes_load(const int n, Newtoncotes **p);
  void newtoncotes_check(const Newtoncotes *p);
  void newtoncotes_free(Newtoncotes *p);
#ifdef __cplusplus
}
#endif

#endif /* NEWTONCOTES_H */
