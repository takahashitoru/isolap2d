#ifndef BINCOEF_H
#define BINCOEF_H

#include <stdio.h>
#include <stdlib.h>

#ifndef BINCOEF
#define BINCOEF(m, b, n, k) ( (b)[(m) * (n) + (k)] )
#else
#error Fail to define BINCOEF
#endif

#ifdef __cplusplus
extern "C" {
#endif
  double *bincoef_malloc(const int m);
  void bincoef_free(double *b);
  void bincoef(const int m, double *b);
#ifdef __cplusplus
}
#endif


#endif /* BINCOEF_H */
