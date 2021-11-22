#ifndef ISOLAP2D_H
#define ISOLAP2D_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "envs.h"
#include "opts.h"

#include "nbie.h"
//140225#include "hnbie.h"
//140225#include "cbie.h"
#include "double2.h"
#include "mathmacros.h"
#include "debugmacros.h"
#include "timer.h"
#include "dgmres.h" // GMRES
#include "gaussone.h"
#include "isolap2d_defs.h"
#include "dxf.h"
#include "options.h"

//140225#define NULL_UINC 0
#include "uinc.h"


/* Set up for FMM */
#include <complex.h>
#include "cell.h"
#include "quadtree.h"
#include "dcomplex.h"
#define CMPLX(a, b) ( (a) + I * (b) )
#include "bincoef.h"

/* etc */
#include "circle.h"
#include "hexa.h"
#include "cshape.h"
#include "sshape.h"


#ifdef __cplusplus
extern "C" {
#endif
  //120124  void isolap2d_input(Isolap2d *p);
  void isolap2d_input_check_cp(Isolap2d *p);
  void isolap2d_input_dxf(char *dxffile, Isolap2d *p);

  void isolap2d_ld(Isolap2d *p);
  void isolap2d_ld_check(const Isolap2d *p);
  void isolap2d_ld_free(Isolap2d *p);


  //140225  void mkuinc(const int ncol, const int ninc, const double2 *xc, double *bvec);
  //140225  void isolap2d_mkuinc(const Isolap2d *p, const int ncol, double *bvec);

  int isolap2d_comp_num_collocation_points(const int ncc, const Closedcurve *cc);
  void isolap2d_comp_collocation_points(const int ncc, const Closedcurve *cc, double *tc, double2 *xc, int *ccc, int *xoff);
  void isolap2d_check_collocation_points(const int nnn, const double2 *xc);


  void isolap2d_direct(const double eps, const GaussOne *g, const int p, const int n, const double *t, const double2 *cp, const double *tc, const double2 *xc, const int icc, const int jcc, const int ix, const int j, double *dval);
  void isolap2d_matvec(const int N, const void *A, const double *xvec, double *yvec);  


  void isolap2d_comp_freeterm(Isolap2d *p, const int ncol);
  int isolap2d_comp_num_basis_center(const int ncc, const Closedcurve *cc);
  //140318  void isolap2d_comp_basis_center(const int ncc, const Closedcurve *cc, double2 *xb, int *bcc, int *bj);
  //140318  void isolap2d_check_basis_center(const int nbase, const double2 *xb);
  void isolap2d_comp_basis_center(Isolap2d *p);
  void isolap2d_free_basis_center(Isolap2d *p);
  void isolap2d_check_basis_center(const Isolap2d *p);

  void isolap2d_msolve(const int n, const void *M, const double *r, double *z);

  void isolap2d_output(const int ncc, const Closedcurve *cc, const double *tc, const double2 *xc, const double *xvec);
  void isolap2d_error(const int ncc, const Closedcurve *cc, const double *tc, const double2 *xc, const double *xvec);

  void isolap2d_setup_root(const int nbase, const int ncol, const double ratio, const double2 *xb, const double2 *xc, double2 *rotcnt, double *rotlen);
  //120124  void isolap2d_make_tree(Isolap2d *p, const double2 rotcnt, const double rotlen, const int nbase, const int ncol);
  void isolap2d_make_tree(Isolap2d *p, const double2 rotcnt, const double rotlen,
			  const int nbasis, const int ncol, const double2 *xb, const double2 *xc);
  void isolap2d_free_tree(Isolap2d *p);
  void isolap2d_chtree(Isolap2d *p);

  void isolap2d_mknf(Isolap2d *p);
  void isolap2d_mknf_check(const Isolap2d *p);
  void isolap2d_mknf_free(Isolap2d *p);

  void isolap2d_matvec_fmm(const int N, const void *A, const double *xvec, double *yvec);
  void isolap2d_incal_fmm(Isolap2d *p, const int nin, const double2 *xin, const double *xvec, double *yvec);

  void isolap2d_incal(const Isolap2d *p, const int nin, const double2 *xin, const double *xvec, double *yvec);
  void isolap2d_incal_input(const Isolap2d *p, const char *infile, int *nin, double2 **xin);
  void isolap2d_incal_input_check(const int nin, const double2 *xin);
  void isolap2d_incal_output(const int nin, const double2 *xin, const double *yvec);

  int isolap2d_comp_num_element_center(const int ncc, const Closedcurve *cc);
  void isolap2d_comp_element_center(Isolap2d *p);
  void isolap2d_free_element_center(Isolap2d *p);
  void isolap2d_check_element_center(const Isolap2d *p);

#ifdef __cplusplus
}
#endif

#endif /* ISOLAP2D_H */
