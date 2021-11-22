#ifndef DGMRES_H
#define DGMRES_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//120116#ifndef TRUE
//120116#define TRUE (1)
//120116#endif
//120116
//120116#ifndef FALSE
//120116#define FALSE (0)
//120116#endif

#ifndef DGMRES_TRUE
#define DGMRES_TRUE (1)
#else
#error DGMRES_TRUE is already defined.
#endif

#ifndef DGMRES_FALSE
#define DGMRES_FALSE (0)
#else
#error DGMRES_FALSE is already defined
#endif

#ifndef DGMRES_SUCCESS
#define DGMRES_SUCCESS (0)
#else
#error DGMRES_SUCCESS is already defined
#endif

#ifndef DGMRES_FAIL
#define DGMRES_FAIL (1)
#else
#error DGMRES_FAIL is already defined
#endif

#ifndef DGMRES_KEEP_RESIDUALS
#define DGMRES_KEEP_RESIDUALS 3
#endif

#ifndef DGMRES_DUMMY_RESIDUAL
#define DGMRES_DUMMY_RESIDUAL (- 1.0) // smaller than zero
#else
#error DGMRES_DUMMY_RESIDUAL is already defined.
#endif


#ifdef __cplusplus
extern "C" {
#endif
  void dgmres_info();
  int dgmres(const int n, const void *A, const void *M, double *x, const double *b,  
	     const int m, const double tol, const int maxiter, const int nullx0,
	     void (*matvec)(const int, const void *, const double *, double *),  
	     void (*msolve)(const int, const void *, const double *, double *),  
	     double *resid, int *niter, int *nmv, int *nms);
#ifdef __cplusplus
}
#endif

#endif /* DGMRES_H */
