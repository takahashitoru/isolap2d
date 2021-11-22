#include "isolap2d.h"

/*
  Usage for check
  1) make with with MYDEBUG and CIRCLE flag
  2) ./a.out > gomi
  3) gnuplot
  % plot'gomi'u 1:2 w lp,''u 1:3 w lp,''u 1:4 w lp,''u 1:5 w lp,''u 1:6 w lp,''u 1:7 w lp,''u 1:8 w lp,''u 1:9 w lp  <-- draw the bases when n=8
  4) gnuplot
  % splot 'boundary.txt' w l, 'knot.txt' w p, 'cp.txt' w p, 'collocation.txt' w p, 'basis.txt' w p
  % unset label; load 'knot.label'; load 'cp.label'; load 'collocation.label'; load 'basis.label'
*/

int main (int argc, char **argv)
{
  /* Option */
  char *dxffile, *infile;
  isolap2d_options(argc, argv, &dxffile, &infile);

  /* Check */    
  envs();
  opts();
  
  /* Set timer */
  timerType *timer_all;
  allocTimer(&timer_all);
  initTimer(timer_all);
  startTimer(timer_all);

  /* Allocate the system pointer */
  Isolap2d *p = (Isolap2d *)malloc(sizeof(Isolap2d));
  
  /* Define n, p, and cp */
  if (dxffile != NULL) { // if dxffile is given by the option "-d"
    isolap2d_input_dxf(dxffile, p);
  } else {
#if defined(HEXA)
    isolap2d_input_hexa_info();
    isolap2d_input_hexa(p); // hexagonal(s)
#elif defined(CSHAPE)
    isolap2d_input_cshape_info();
    isolap2d_input_cshape(p); // c-shape(s)
#elif defined(SSHAPE)
    isolap2d_input_sshape_info();
    isolap2d_input_sshape(p); // s-shape(s)
#else
    isolap2d_input_circle_info();
    isolap2d_input_circle(p); // circle(s)
#endif
  }
#if defined(MYDEBUG)
    isolap2d_input_check_cp(p);
#endif

  /* Load parameters */
  isolap2d_ld(p);
  isolap2d_ld_check(p);
  
  /* Load Gaussian quadrature formula */
  GaussOne_load(p->params->ngauss, &(p->gauss));
#if defined(MYDEBUG)
  GaussOne_check(p->gauss);
#endif

  /* Compute knots */
  for (int icc = 0; icc < p->ncc; icc ++) {
    (p->cc)[icc].t = closedcurve_knot_alloc((p->cc)[icc].p, (p->cc)[icc].n); // t[0:n+p]
    closedcurve_knot_comp((p->cc)[icc].p, (p->cc)[icc].n, (p->cc)[icc].t);
  }
#if defined(MYDEBUG)
  closedcurve_knot_check(p->ncc, p->cc);
  closedcurve_draw_boundary(p->ncc, p->cc);
#endif

  /* Compute collocation points */
  const int ncol = isolap2d_comp_num_collocation_points(p->ncc, p->cc);
  INFO("ncol=%d\n", ncol);
  p->tc = (double *)malloc(ncol * sizeof(double)); // tc[0:ncol-1]
  p->xc = (double2 *)malloc(ncol * sizeof(double2)); // xc[0:ncol-1]
  p->ccc = (int *)malloc(ncol * sizeof(int)); // ccc[0:ncol-1]
  p->xoff = (int *)malloc(p->ncc * sizeof(int)); // xoff[0:ncc-1]
  isolap2d_comp_collocation_points(p->ncc, p->cc, p->tc, p->xc, p->ccc, p->xoff);
#if defined(MYDEBUG)
  isolap2d_check_collocation_points(ncol, p->xc);
#endif

#if defined(FMM)
  //140318  /* Compute centres of basis functions */
  //140318  const int nbasis = isolap2d_comp_num_basis_center(p->ncc, p->cc);
  //140318  INFO("nbasis=%d\n", nbasis);
  //140318  p->xb = (double2 *)malloc(nbasis * sizeof(double2)); // xb[0:nbasis-1]
  //140318  p->bcc = (int *)malloc(nbasis * sizeof(int)); // bcc[0:nbasis-1]
  //140318  p->bj = (int *)malloc(nbasis * sizeof(int)); // bj[0:nbasis-1]
  //140318 isolap2d_comp_basis_center(p->ncc, p->cc, p->xb, p->bcc, p->bj);
  /* Compute centres of bases or elements */
#if defined(ORIGINAL)
  isolap2d_comp_basis_center(p);
#else
  isolap2d_comp_element_center(p);
#endif
#if defined(MYDEBUG)
  //140318  isolap2d_check_basis_center(nbasis, p->xb);
#if defined(ORIGINAL)
  isolap2d_check_basis_center(p);
#else
  isolap2d_check_element_center(p);
#endif
#endif
  
  /* Build the quadtree */
  double2 rotcnt;
  double rotlen;
  //140318  isolap2d_setup_root(nbasis, ncol, p->params->ratio, p->xb, p->xc, &rotcnt, &rotlen);
  //140318  isolap2d_make_tree(p, rotcnt, rotlen, nbasis, ncol, p->xb, p->xc);
#if defined(ORIGINAL)
  const int nsource = isolap2d_comp_num_basis_center(p->ncc, p->cc);
#else
  const int nsource = isolap2d_comp_num_element_center(p->ncc, p->cc);
#endif
  isolap2d_setup_root(nsource, ncol, p->params->ratio, p->xb, p->xc, &rotcnt, &rotlen);
  isolap2d_make_tree(p, rotcnt, rotlen, nsource, ncol, p->xb, p->xc);
#if defined(MYDEBUG)
  isolap2d_chtree(p);
#endif

  /* Compute neighbour and interaction lists */
  isolap2d_mknf(p);
#if defined(MYDEBUG)
  isolap2d_mknf_check(p);
#endif

#endif /* FMM */

  /* Compute free term */
  p->freeterm = (double *)malloc(ncol * sizeof(double)); // freeterm[0:ncol-1];
#if defined(ORIGINAL) || defined(CBIE)
  isolap2d_comp_freeterm(p, ncol);
#endif

  /* Compute b-vector */
  double *bvec = (double *)calloc(ncol, sizeof(double)); // bvec[0:ncol-1]; initialise
#if defined(ORIGINAL) || defined(CBIE)
  isolap2d_mkuinc(p, ncol, bvec); // contribution from the incidental field unless ninc=0
#else
  isolap2d_mkqinc(p, ncol, bvec); // contribution from the incidental field unless ninc=0
#endif

  /* Solve Ax=b by GMRES */
  double *xvec = (double *)malloc(ncol * sizeof(double)); // xvec[0:ncol-1]
  for (int i = 0; i < ncol; i ++) {
    xvec[i] = 0.0; // Initial guess is set to zero (then you can set initx to 1)
  }

  dgmres_info(); // Print macros defined in dgmres
  
  int niter, nmv, nms; // #s of iterations, matvec, and msolve until convergence
  double resid; // residual
  Params *pa = p->params; // alias

#if defined(FMM)
  int rval = dgmres(ncol, (void *)p, (void *)p, xvec, bvec, pa->maxl, pa->tol, pa->itmax, pa->initx,
		    &isolap2d_matvec_fmm, &isolap2d_msolve, &resid, &niter, &nmv, &nms);
#else
  //140318#if defined(ORIGINAL)
  int rval = dgmres(ncol, (void *)p, (void *)p, xvec, bvec, pa->maxl, pa->tol, pa->itmax, pa->initx,
		    &isolap2d_matvec, &isolap2d_msolve, &resid, &niter, &nmv, &nms);
  //140318#else
  //140318  int rval = dgmres(ncol, (void *)p, (void *)p, xvec, bvec, pa->maxl, pa->tol, pa->itmax, pa->initx,
  //140318  		    &isolap2d_matvec_nbie, &isolap2d_msolve, &resid, &niter, &nmv, &nms);
  //140318#endif
#endif
  if (rval == DGMRES_FAIL) {
    INFO("dgmres failed to converge; niter=%d nmv=%d nms=%d resid=%e\n", niter, nmv, nms, resid);
  }

#if defined(FMM)
  /* Free quadtree memories associated with collocation points */
  isolap2d_mknf_free(p);
  isolap2d_free_tree(p);
#endif

  /* Output the solution */
  isolap2d_output(p->ncc, p->cc, p->tc, p->xc, xvec);

  /* Calculation for internal points */
  if (infile != NULL) {
    int nin;
    double2 *xin;
    isolap2d_incal_input(p, infile, &nin, &xin);
#if defined(MYDEBUG)
    isolap2d_incal_input_check(nin, xin);
#endif

    double *yvec = (double *)malloc(nin * sizeof(double)); // store potentials at internal points

#if defined(FMM)
    double2 rotcnt;
    double rotlen;
    //140318    isolap2d_setup_root(nbasis, nin, p->params->ratio, p->xb, xin, &rotcnt, &rotlen); // ncol => nin, xc => xin
    //140318    isolap2d_make_tree(p, rotcnt, rotlen, nbasis, nin, p->xb, xin); // ncol => nin, xc => xin
    isolap2d_setup_root(nsource, nin, p->params->ratio, p->xb, xin, &rotcnt, &rotlen); // ncol => nin, xc => xin
    isolap2d_make_tree(p, rotcnt, rotlen, nsource, nin, p->xb, xin); // ncol => nin, xc => xin
#if defined(MYDEBUG)
    isolap2d_chtree(p);
#endif
    isolap2d_mknf(p);
#if defined(MYDEBUG)
    isolap2d_mknf_check(p);
#endif
    isolap2d_incal_fmm(p, nin, xin, xvec, yvec);
    isolap2d_mknf_free(p);
    isolap2d_free_tree(p);
#else // conventional
    isolap2d_incal(p, nin, xin, xvec, yvec);
#endif

    mkuinc(nin, p->params->ninc, xin, yvec); // add the external field
    isolap2d_incal_output(nin, xin, yvec);

    free(yvec);
    free(xin);
  }

  /* Stop timer (exclude computation of error) */
  stopTimer(timer_all);
  printTimer(stderr, "timer_all", timer_all);
  freeTimer(&timer_all);

  /* Compute L2-error */
#if defined(CIRCLE)
  if (pa->ninc == 300) {
    isolap2d_error(p->ncc, p->cc, p->tc, p->xc, xvec);
  } else {
    MESG("Unexpected uinc.\n");
  }
#endif


  

  /* Finalise */
  free(xvec);
  free(bvec);
#if defined(FMM)
  //140318  free(p->xb);
  //140318  free(p->bcc);
  //140318  free(p->bj);
#if defined(ORIGINAL)
  isolap2d_free_basis_center(p);
#else
  isolap2d_free_element_center(p);
#endif
#endif
  free(p->freeterm);

  free(p->tc);
  free(p->xc);
  free(p->ccc);
  free(p->xoff);
  GaussOne_free(p->gauss);
  isolap2d_ld_free(p);
  for (int icc = 0; icc < p->ncc; icc ++) {
    closedcurve_knot_free((p->cc)[icc].t);
    closedcurve_cp_free((p->cc)[icc].cp);
  }
  closedcurve_free(p->cc);

  free(p);

  return(EXIT_SUCCESS);
}
