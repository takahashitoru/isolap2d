#ifndef DOUBLE2_H
#define DOUBLE2_H

#include <math.h>

typedef struct {
  double x;
  double y;
} double2;

#ifdef __cplusplus
extern "C" {
#endif
  double dot2(const double2 u, const double2 v);
  double norm2(const double2 u);  
  double2 add2(const double2 u, const double2 v);
  double2 sub2(const double2 u, const double2 v);
  double2 scale2(const double a, const double2 u);  
  double vdot2(const double2 u, const double2 v);
  double2 set2(const double x, const double y);
#ifdef __cplusplus
}
#endif

#endif /* DOUBLE2_H */
