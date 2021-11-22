#ifndef NBIE_H
#define NBIE_H

#include "double2.h"
#include "gaussone.h"
#include "isolap2d_defs.h"

#include "uinc.h"

#define SMALL_NUMBER ( 1.0e-3 )
//#define SMALL_NUMBER ( 0.0 )

#ifdef __cplusplus
extern "C" {
#endif
  void isolap2d_mkqinc(const Isolap2d *p, const int ncol, double *bvec);  
  //140318  void isolap2d_matvec_nbie(const int N, const void *A, const double *xvec, double *yvec);
  void matvec_nbie(const GaussOne *g, const int ncc, const Closedcurve *cc, const int ncol, const double *tc, const double2 *xc, const double *freeterm, const double *xvec, double *yvec);
  
  //140318  double indefinite_integral_IB(double a, const double b, const double c, const double x, const double D);
  //140318  double integrand1_1(const double a, const double b, const double c, const double IB, const double x); // from formula
  //140318  double integrand2_1(const double a, const double b, const double c, const double IB, const double x); // from formula
  //140318  double integrand0_2(const double a, const double b, const double c, const double IB, const double x, const double D); // from formula
  //140318  double integrand1_2(const double a, const double b, const double c, const double IB, const double x, const double D); // from formula
  //140318  double integrand2_2(const double a, const double b, const double c, const double IB, const double x, const double D); // from formula

  //140325  void isolap2d_direct_nbie_fmm(const GaussOne *g, const int p, const int n, const double *t, const double2 *cp, const double *tc, const double2 *xc, const int icc, const int jcc, double *dval, const double *xvec, const int ioff, const int joff, const int i, const int j, const int pi, const int ni, const double *ti, const double2 *cpi);

  void isolap2d_direct_nbie(const GaussOne *g, const int p, const int n, const double *t, const double2 *cp, const int icc, const int jcc, double *dval, const double *xvec, const int ioff, const int joff, const int i, const double2 x, const double2 nx, const int j);

#ifdef __cplusplus
}
#endif

#endif /* NBIE_H */
