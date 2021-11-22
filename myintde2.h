#ifndef MYINTDE2_H
#define MYINTDE2_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef MYINTDE_LENAW
#define MYINTDE_LENAW ( 8000 ) // this value is from intde2t.c
#endif
#ifndef MYINTDE_TINY
#define MYINTDE_TINY ( 1.0e-307 )  // this value is from intde2t.c
#endif
#ifndef MYINTDE_EPS
#define MYINTDE_EPS ( 1.0e-15 )  // this value is from intde2t.c
#endif


#if(1) // Define your own structure

#include "double2.h"

typedef struct {
  int p;
  int n;
  double *t;
  double2 *cp;
  double2 x;
  double2 nx;
  double *xvec;
  int joff;
  int k;
} IntdeParams;

#else // This is an example for test_myintde2.c

typedef struct {
  double a;
  double b;
  double c;
} IntdeParams;

#endif

typedef struct {
  int lenaw; // length of aw (int)
  double *aw;
  double tiny; // minimum value that 1/tiny does not overflow (double)
  double eps; // relative error requested (double)
  IntdeParams *params;
} Intde;

#ifdef __cplusplus
extern "C" {
#endif
  Intde *myintde_malloc();
  void myintde_init(Intde *p);
  void myintde(Intde *p, double (*f)(double, Intde *), double a, double b, double *i, double *err);
  void myintde_free(Intde *p);
#ifdef __cplusplus
}
#endif


#endif /* MYINTDE2_H */
