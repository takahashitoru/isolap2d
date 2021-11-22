#include "isolap2d.h"


#if !defined(FUKUI)
static void comp_fi(const int m, const dcomplex z, dcomplex *fi)
{
  fi[0] = 1.0;
  for (int n = 1; n < m; n ++) {
    fi[n] = fi[n - 1] * z / n;
  }
}
#endif


static void y2m(const GaussOne *g, const Closedcurve *cc, const int *barray, const int *bcc, const int *bj, const int *xoff,
		Cell *c, const int nterm, const int icell, const double *xvec)
{
  //140318  /* Constant */
  //140318  const double pi2i = 1.0 / (2.0 * PI);

  /* Set up for Guassian quadrature */
  const int ng = g->n;
  const double *gx = g->x;
  const double *gw = g->w;

  /* Loop over bases in this cell */
  for (int kb = c[icell].bhead; kb <= c[icell].btail; kb ++) {
    
    /* Obtain the index of this basis, i.e., b_{j,p} */
    const int j = bj[barray[kb]]; // global to local
    const int jcc = bcc[barray[kb]]; // index of closed curve
    const int p = cc[jcc].p; // degree
    const int n = cc[jcc].n; // number of control points
    const double *t = cc[jcc].t; // knot values
    const double2 *cp = cc[jcc].cp; // control points
    const int joff = xoff[jcc];
      
    /* Loop over the intervals in the support of b_{j,p} */
    for (int k = 0; k <= p; k ++) {
      
      /* Compute t such as on the boundary, i.e., t in [t_p,t_n] */
      if (j + k >= p && j + k + 1 <= n) {
	
	/* Compute the Jacobian wrt the Gaussian quadrature (this is
	   not the Jacobian wrt the curve parameter) */
	const double Jg = (t[j + k + 1] - t[j + k]) / 2.0;
      
	/* Obtain the potential u_j by considering the wrapping */
	const double uj = xvec[joff + j % (n - p)];

	///////////////////////////////////////////////////////////
	//	DBG("icell=%d kb=%d j=%d k=%d uj=%e\n", icell, kb, j, k, uj);
	///////////////////////////////////////////////////////////

	/* Loop over Gaussian abscissa */
	for (int ig = 0; ig < ng; ig ++) {
	
	  /* Parameter corresponding to the ig th abscissa */
	  const double tnow = ((1.0 - gx[ig]) * t[j + k] + (1.0 + gx[ig]) * t[j + k + 1]) / 2.0;
	
	  /* Compute the point y and its derivative dy/dt */
	  const double2 y = closedcurve_position(p, n, t, cp, tnow);
	  const double2 dy = closedcurve_position_derivative(p, n, t, cp, tnow);
	  const dcomplex T = CMPLX(dy.x, dy.y);

	  /* Compute the vector from centre to y */
	  const dcomplex z = CMPLX(y.x - c[icell].center.x, y.y - c[icell].center.y);

	  /* Compute b_{j,p}(tnow) */
	  const double bj = closedcurve_basis(p, n, t, j, tnow);

#if defined(FUKUI)
	  /* Compute moment */
	  c[icell].mcoe[0] = 0.0; // always zero
	  dcomplex ztmp = 1.0; // z^{m-1}
	  for (int m = 1; m < nterm; m ++) {
	    //140318	    c[icell].mcoe[m] += pi2i * gw[ig] * (- I) * T * ztmp * bj * Jg * uj;
	    c[icell].mcoe[m] += PI2I * gw[ig] * (- I) * T * ztmp * bj * Jg * uj;
	    ztmp *= z;
	  }
#else
	  /* Compute I function */
	  dcomplex *fi = (dcomplex *)malloc((nterm - 1) * sizeof(dcomplex)); // fi[0:nterm-2]
	  comp_fi(nterm - 1, z, fi);
	  
	  /* Compute moment */
	  c[icell].mcoe[0] = 0.0; // always zero
	  for (int m = 1; m < nterm; m ++) {
	    //140318	    c[icell].mcoe[m] += pi2i * gw[ig] * (- I) * T * fi[m - 1] * bj * Jg * uj;
	    c[icell].mcoe[m] += PI2I * gw[ig] * (- I) * T * fi[m - 1] * bj * Jg * uj;
	  }
	  
	  free(fi);
#endif

	} // ig
      }
    } // k
  } // b

}



static void y2m_nbie(const GaussOne *g, const Closedcurve *cc, const int *barray, const int *bcc, const int *bj, const int *xoff, Cell *c, const int nterm, const int icell, const double *xvec)
{
  /* Set up for Guassian quadrature */
  const int ng = g->n;
  const double *gx = g->x;
  const double *gw = g->w;

  /* Loop over isogeometric elements in this cell */
  for (int kb = c[icell].bhead; kb <= c[icell].btail; kb ++) {
    
    /* Perform the integral over this element, i.e., I_j:=[t_{j+p}, t_{j+p+1}) */
    const int j = bj[barray[kb]]; // from global index of this element to local one
    const int jcc = bcc[barray[kb]]; // index of closed curve
    const int p = cc[jcc].p; // degree
    const int n = cc[jcc].n; // number of control points
    const double *t = cc[jcc].t; // knot values
    const double2 *cp = cc[jcc].cp; // control points
    const int joff = xoff[jcc];
      
    /* Compute the Jacobian wrt the Gaussian quadrature (this is not
       the Jacobian wrt the curve parameter) */
    const double Jg = (t[j + p + 1] - t[j + p]) / 2.0;
    
    /* Loop over Gaussian abscissa */
    for (int ig = 0; ig < ng; ig ++) {
	
      /* Parameter corresponding to the ig th abscissa */
      const double tnow = ((1.0 - gx[ig]) * t[j + p] + (1.0 + gx[ig]) * t[j + p + 1]) / 2.0;
	
      /* Compute the point y and its derivative dy/dt */
      const double2 y = closedcurve_position(p, n, t, cp, tnow);
      const double2 dy = closedcurve_position_derivative(p, n, t, cp, tnow);
      const dcomplex T = CMPLX(dy.x, dy.y);

      /* Compute the vector from centre to y */
      const dcomplex z = CMPLX(y.x - c[icell].center.x, y.y - c[icell].center.y);

      /* Compute u as the sum of (p+1) bases, viz. b_{j},...,b_{j+p} */
      double uval = 0.0;
      for (int k = 0; k <= p; k ++) {
	const double bsp = closedcurve_basis(p, n, t, j + k, tnow);
	uval += bsp * xvec[joff + (j + k) % (n - p)];
      }

#if defined(FUKUI)
      /* Compute moment */
      c[icell].mcoe[0] = 0.0; // always zero
      dcomplex ztmp = 1.0; // z^{m-1}
      for (int m = 1; m < nterm; m ++) {
	c[icell].mcoe[m] += PI2I * gw[ig] * (- I) * T * ztmp * Jg * uval;
	ztmp *= z;
      }
#else
      /* Compute I function */
      dcomplex *fi = (dcomplex *)malloc((nterm - 1) * sizeof(dcomplex)); // fi[0:nterm-2]
      comp_fi(nterm - 1, z, fi);
	  
      /* Compute moment */
      c[icell].mcoe[0] = 0.0; // always zero
      for (int m = 1; m < nterm; m ++) {
	c[icell].mcoe[m] += PI2I * gw[ig] * (- I) * T * fi[m - 1] * Jg * uval;
      }
	  
      free(fi);
#endif

    } // ig

  } // kb
}



static void m2m(Cell *c, const int nterm, const int icell)
{
#if defined(FUKUI)
  double *bi = bincoef_malloc(nterm - 1); // bi[0:nterm-2][0:nterm-2]
  bincoef(nterm - 1, bi);
#endif

  for (int ic = 0; ic < 4; ic ++) {
    const int jcell = c[icell].child[ic];

    if (jcell != NULL_CELL) {

      const dcomplex z = CMPLX(c[jcell].center.x - c[icell].center.x, c[jcell].center.y - c[icell].center.y);

#if defined(FUKUI)
      dcomplex *zpow = (dcomplex *)malloc(nterm * sizeof(dcomplex));
      dcomplex ztmp = 1.0;
      for (int n = 0; n < nterm; n ++) {
	zpow[n] = ztmp;
	ztmp *= z;
      }
      c[icell].mcoe[0] += c[jcell].mcoe[0];
      for (int n = 1; n < nterm; n ++) {
	c[icell].mcoe[n] += zpow[n] / n * c[jcell].mcoe[0];
	for (int m = 1; m <= n; m ++) {
	  c[icell].mcoe[n] += BINCOEF(nterm - 1, bi, n - 1, m - 1) * zpow[n - m] * c[jcell].mcoe[m];
	}
      }
      free(zpow);
#else
      dcomplex *fi = (dcomplex *)malloc(nterm * sizeof(dcomplex)); // fi[0:nterm-1]
      comp_fi(nterm, z, fi);
      for(int n = 0; n < nterm; n ++) {
	dcomplex ztmp = 0.0;
	for (int k = 0; k <= n; k ++) {
	  ztmp += c[jcell].mcoe[k] * fi[n - k];
	}
	c[icell].mcoe[n] += ztmp;
      }
      free(fi);
#endif

    }
  }

#if defined(FUKUI)
  bincoef_free(bi);
#endif
}


#if !defined(FUKUI)
static void comp_fo(const int m, const dcomplex z, dcomplex *fo)
{
  //////////////////////
  ASSERT(cabs(z) > 0.0);
  ASSERT(m > 1);
  //////////////////////
  fo[0] = - clog(z);
  fo[1] = 1.0 / z;
  for (int n = 2; n < m; n ++) {
    fo[n] = fo[n - 1] * (n - 1) / z;
  }
}
#endif

static void m2l(Cell *c, const int nterm, const int icell)
{
#if defined(FUKUI)
  double *bi = bincoef_malloc(2 * nterm - 2); // bi[0:2*nterm-3][0:2*nterm-3]
  bincoef(2 * nterm - 2, bi);
#endif

  for (int iinter = 0; iinter < c[icell].ninter; iinter ++) {

    const int jcell = c[icell].interactions[iinter];

    const dcomplex z = CMPLX(c[icell].center.x - c[jcell].center.x, c[icell].center.y - c[jcell].center.y);

#if defined(FUKUI)
    const dcomplex zinv = 1.0 / z;
    dcomplex zsum = clog(zinv) * c[jcell].mcoe[0];
    dcomplex ztmp = zinv; // 1/z^m
    for (int m = 1; m < nterm; m ++) {
      zsum += ztmp * c[jcell].mcoe[m];
      ztmp *= zinv;
    }
    c[icell].lcoe[0] += zsum;
    ztmp = - zinv; // (-1)^n/z^n
    for (int n = 1; n < nterm; n ++) {
      zsum = c[jcell].mcoe[0] / n;
      dcomplex ztmp2 = zinv; // 1/z^m
      for (int m = 1; m < nterm; m ++) {
	zsum += BINCOEF(2 * nterm - 2, bi, n + m - 1, m - 1) * ztmp2 * c[jcell].mcoe[m];
	ztmp2 *= zinv;
      }
      c[icell].lcoe[n] += ztmp * zsum;
      ztmp *= - zinv;
    }
#else
    dcomplex *fo = (dcomplex *)malloc((2 * nterm - 1) * sizeof(dcomplex)); // fo[0:2*nterm-2]  
    comp_fo(2 * nterm - 1, z, fo);
    for (int n = 0; n < nterm; n ++) {
      for (int k = 0; k < nterm; k ++) {
	c[icell].lcoe[n] += c[jcell].mcoe[k] * fo[n + k];
      }
    }
    free(fo);
#endif

  } // iinter

#if defined(FUKUI)
  bincoef_free(bi);
#endif

}


static void l2l(Cell *c, const int nterm, const int icell)
{
#if defined(FUKUI)
  double *bi = bincoef_malloc(nterm); // bi[0:nterm-1][0:nterm-1]
  bincoef(nterm, bi);
#endif

  for (int i = 0; i < 4; i ++) {
    const int jcell = c[icell].child[i];

    if (jcell != NULL_CELL) {

#if defined(FUKUI)
      const dcomplex z = CMPLX(c[jcell].center.x - c[icell].center.x, c[jcell].center.y - c[icell].center.y);
      dcomplex *zpow = (dcomplex *)malloc(nterm * sizeof(dcomplex)); // zpow[0:nterm-1]
      dcomplex ztmp = 1.0; // z^n
      for (int n = 0; n < nterm; n ++) {
	zpow[n] = ztmp;
	ztmp *= z;
      }
      for (int n = 0; n < nterm; n ++) {
	dcomplex zsum = 0.0;
	for (int m = n; m < nterm; m ++) {
	  zsum += BINCOEF(nterm, bi, m, n) * zpow[m - n] * c[icell].lcoe[m];
	}
	c[jcell].lcoe[n] += zsum;
      }
      free(zpow);
#else
      const dcomplex z = CMPLX(c[icell].center.x - c[jcell].center.x, c[icell].center.y - c[jcell].center.y);
      dcomplex *fi = (dcomplex *)malloc(nterm * sizeof(dcomplex)); // fi[0:nterm-1]
      comp_fi(nterm, z, fi);
      for (int n = 0; n < nterm; n ++) {
	dcomplex zsum = 0.0;
	for (int k = n; k < nterm; k ++) {
	  zsum += c[icell].lcoe[k] * fi[k - n];
	}
	c[jcell].lcoe[n] += zsum;
      }
      free(fi);
#endif

    }
  }

#if defined(FUKUI)
  bincoef_free(bi);
#endif

}


static void l2x(const int *carray, const double2 *xc, Cell *c, const int nterm, const int icell, double *yvec)
{
  for (int kx = c[icell].chead; kx <= c[icell].ctail; kx ++) {

    const int ix = carray[kx]; // obtain a collocation point in this icell

#if defined(FUKUI)
    const dcomplex z = CMPLX(xc[ix].x - c[icell].center.x, xc[ix].y - c[icell].center.y);
    dcomplex ztmp = 1.0; // z^l
    for (int l = 0; l < nterm; l ++) {
      yvec[ix] += creal(c[icell].lcoe[l] * ztmp);
      ztmp *= z;
    }
#else
    const dcomplex z = CMPLX(c[icell].center.x - xc[ix].x, c[icell].center.y - xc[ix].y);
    dcomplex *fi = (dcomplex *)malloc(nterm * sizeof(dcomplex)); // fi[0:nterm-1]
    comp_fi(nterm, z, fi);
    for (int l = 0; l < nterm; l ++) {
      yvec[ix] += creal(c[icell].lcoe[l] * fi[l]);
    }
    free(fi);
#endif
  }
}


static void l2x_nbie(const int *carray, const Closedcurve *cc, const int *ccc, const double *tc, const double2 *xc, Cell *c, const int nterm, const int icell, double *yvec)
{
  const dcomplex *lcoe = c[icell].lcoe;

  for (int kx = c[icell].chead; kx <= c[icell].ctail; kx ++) {

    const int ix = carray[kx]; // obtain a collocation point in this icell
    const int icc = ccc[ix]; // closed curve
    const int pi = cc[icc].p; // degree
    const int ni = cc[icc].n; // number of control points
    const double *ti = cc[icc].t; // knot values
    const double2 *cpi = cc[icc].cp; // control points

    const double2 dx = closedcurve_position_derivative(pi, ni, ti, cpi, tc[ix]);
    double2 nx; // outward normal vector at the collocation point
    nx.x = dx.y / norm2(dx);
    nx.y = - dx.x / norm2(dx);

    dcomplex zsum = 0.0;

#if defined(FUKUI)
    const dcomplex z = CMPLX(xc[ix].x - c[icell].center.x, xc[ix].y - c[icell].center.y);
    dcomplex ztmp = 1.0; // z^{l-1} where l=1
    for (int l = 1; l < nterm; l ++) {
      zsum += lcoe[l] * ztmp * l; // L_l * z^{l-1} * l
      ztmp *= z;
    }
#else
    const dcomplex z = CMPLX(c[icell].center.x - xc[ix].x, c[icell].center.y - xc[ix].y);
    dcomplex *fi = (dcomplex *)malloc((nterm - 1) * sizeof(dcomplex)); // fi[0:nterm-2]
    comp_fi(nterm - 1, z, fi);
    for (int l = 1; l < nterm; l ++) {
      zsum += lcoe[l] * (- fi[l - 1]);
    }
    free(fi);
#endif

    yvec[ix] += nx.x * creal(zsum) - nx.y * cimag(zsum);

  }
}


static void upward(const int nterm, const int maxlev, const int minlev,
		   const int *barray, const Closedcurve *cc, const int *bcc, const int *bj, const int *xoff,
		   const int *levsta, const int *levend, Cell *c, 
		   const GaussOne *g, const double *xvec)
{
  for (int level = maxlev; level >= minlev; level --) {
    
#if defined(_OPENMP)
#pragma omp parallel for
#endif
    for (int icell = levsta[level]; icell <= levend[level]; icell ++) {
      
      if (imleaf(c[icell]) == QUADTREE_TRUE) {

#if defined(ORIGINAL)	
	y2m(g, cc, barray, bcc, bj, xoff, c, nterm, icell, xvec); // basis-wise
#else
	y2m_nbie(g, cc, barray, bcc, bj, xoff, c, nterm, icell, xvec); // element-wise
#endif

      } else {

	m2m(c, nterm, icell); // from children

      }

    }
  }

}  


static void downward(const int nterm, const int maxlev, const int minlev,
		     const int *barray, const int *carray, const int *levsta, const int *levend, Cell *c,
		     const double eps, const GaussOne *g,
		     const Closedcurve *cc, const int *bcc, const int *bj, const int *xoff, const int *ccc,
		     const double *tc, const double2 *xc, const double *freeterm, const double *xvec, double *yvec)
{

  for (int level = minlev; level <= maxlev; level ++) {

#if defined(_OPENMP)
#pragma omp parallel for
#endif
    for (int icell = levsta[level]; icell <= levend[level]; icell ++) {
      
      m2l(c, nterm, icell);

      const int iflag = imleaf(c[icell]);
      
      if (iflag == QUADTREE_TRUE) {

#if defined(ORIGINAL) || defined(CBIE)
	l2x(carray, xc, c, nterm, icell, yvec);
#else
	l2x_nbie(carray, cc, ccc, tc, xc, c, nterm, icell, yvec);
#endif

#if defined(ORIGINAL) || defined(CBIE)
	for (int kx = c[icell].chead; kx <= c[icell].ctail; kx ++) {
	  const int ix = carray[kx];
	  const int icc = ccc[ix];
	  const int pi = cc[icc].p;
	  const int ni = cc[icc].n;
	  const double *ti = cc[icc].t;
	  yvec[ix] += freeterm[ix] * closedcurve_potential(pi, ni, ti, &(xvec[xoff[icc]]), tc[ix]); // add freeterm
	}
#endif

      } else {
	
	l2l(c, nterm, icell); // to children
	
      }
      
      for (int ineigh = 0; ineigh < c[icell].nneigh; ineigh ++) {
	const int jcell = c[icell].neighbors[ineigh];
	const int jflag = imleaf(c[jcell]);
	if (iflag == QUADTREE_TRUE || jflag == QUADTREE_TRUE) {
	  for (int kx = c[icell].chead; kx <= c[icell].ctail; kx ++) {
	    const int ix = carray[kx]; // global index of collocation point
	    const int icc = ccc[ix];
#if !defined(ORIGINAL)
	    const int pi = cc[icc].p;
	    const int ni = cc[icc].n;
	    const double *ti = cc[icc].t;
	    const double2 *cpi = cc[icc].cp;
	    const int ioff = xoff[icc];
	    const int i = ix - ioff; // local index of collocation point
	    const double2 x = xc[ix]; // collocation point
	    const double2 dx = closedcurve_position_derivative(pi, ni, ti, cpi, tc[ix]);
	    const double2 nx = scale2(1.0 / norm2(dx), set2(dx.y, - dx.x)); // outward normal vector at x
#endif
	    for (int kb = c[jcell].bhead; kb <= c[jcell].btail; kb ++) {
	      //140325	      const int j = bj[barray[kb]]; // local index of basis or element in jcell
	      //140325	      const int jcc = bcc[barray[kb]];
	      const int iy = barray[kb]; // global index of basis or element
	      const int j = bj[iy]; // local index of basis or element
	      const int jcc = bcc[iy];
	      const int p = cc[jcc].p;
	      const int n = cc[jcc].n;
	      const double *t = cc[jcc].t;
	      const double2 *cp = cc[jcc].cp;
	      const int joff = xoff[jcc];
	      //////////////////////////////////
	      ASSERT(iy == j + joff);
	      //////////////////////////////////
	      double dval;
#if defined(ORIGINAL)
	      isolap2d_direct(eps, g, p, n, t, cp, tc, xc, icc, jcc, ix, j, &dval);
	      yvec[ix] += dval * xvec[joff + j % (n - p)];
#else
	      //140325	      isolap2d_direct_nbie_fmm(g, p, n, t, cp, tc, xc, icc, jcc, &dval, xvec, ioff, joff, i, j, pi, ni, ti, cpi);
	      isolap2d_direct_nbie(g, p, n, t, cp, icc, jcc, &dval, xvec, ioff, joff, i, x, nx, j);
	      yvec[ix] += dval;
#endif
	    }
	  }
	}
      }

      
    } // icell
  } // level

}


static void matvec_fmm(const int ncol, const int nterm, const int maxlev, const int minlev, const int ncell,
		       const int *barray, const int *carray, const int *levsta, const int *levend, Cell *c, 
		       const double eps, const GaussOne *g,
		       const Closedcurve *cc, const int *bcc, const int *bj, const int *xoff, const int *ccc,
		       const double *tc, const double2 *xc, const double *freeterm, const double *xvec, double *yvec)
{
  for (int icell = 0; icell < ncell; icell ++) {
    c[icell].mcoe = (dcomplex *)calloc(nterm, sizeof(dcomplex)); // initialise
    c[icell].lcoe = (dcomplex *)calloc(nterm, sizeof(dcomplex)); // initialise
  }
  
  for(int ix = 0; ix < ncol; ix ++) {
    yvec[ix] = 0; // initialise
  }

  ////////////////////////////////////
  DBG("Enter upward\n");
  ////////////////////////////////////
  upward(nterm, maxlev, minlev,
	 barray, cc, bcc, bj, xoff,
	 levsta, levend, c,
	 g, xvec);
  ////////////////////////////////////
  DBG("Enter downward\n");
  ////////////////////////////////////
  downward(nterm, maxlev, minlev,
	   barray, carray, levsta, levend, c, 
	   eps, g,
	   cc, bcc, bj, xoff, ccc,
	   tc, xc, freeterm, xvec, yvec);

  for (int icell = 0; icell < ncell; icell ++) {
    free(c[icell].mcoe);
    free(c[icell].lcoe);
  }
}


void isolap2d_matvec_fmm(const int N, const void *A, const double *xvec, double *yvec)
{
  Isolap2d *p = (Isolap2d *)A;
  Params *pa = p->params;
  Tree *tree = p->tree;

  matvec_fmm(N, pa->nterm, tree->maxlev, tree->minlev, tree->ncell,
	     tree->barray, tree->carray, tree->levsta, tree->levend, tree->cell, 
	     pa->eps, p->gauss, 
	     p->cc, p->bcc, p->bj, p->xoff, p->ccc,
	     p->tc, p->xc, p->freeterm, xvec, yvec);
}



//static void downward(const int nterm, const int maxlev, const int minlev,
//		     const int *barray, const int *carray, const int *levsta, const int *levend, Cell *c,
//		     const double eps, const GaussOne *g,
//		     const Closedcurve *cc, const int *bcc, const int *bj, const int *xoff, const int *ccc,
//		     const double *tc, const double2 *xc, const double *freeterm, const double *xvec, double *yvec)
static void downward_incal(const int nterm, const int maxlev, const int minlev,
			   const int *barray, const int *carray, const int *levsta, const int *levend, Cell *c,
			   const double eps, const GaussOne *g,
			   const Closedcurve *cc, const int *bcc, const int *bj, const int *xoff,
			   const double *tc, const double2 *xin, const double *xvec, double *yvec)
{

  for (int level = minlev; level <= maxlev; level ++) {

#if defined(_OPENMP)
#pragma omp parallel for // OpenMP DEFINED LOOP WAS PARALLELIZED.
#endif
    for (int icell = levsta[level]; icell <= levend[level]; icell ++) {
      
      m2l(c, nterm, icell);

      const int iflag = imleaf(c[icell]);
      
      if (iflag == QUADTREE_TRUE) {

	//	l2x(carray, xc, c, nterm, icell, yvec);
	l2x(carray, xin, c, nterm, icell, yvec);

	//	for (int kx = c[icell].chead; kx <= c[icell].ctail; kx ++) {
	//	  const int ix = carray[kx]; // collocation point in icell
	//	  const int icc = ccc[ix]; // closed curve
	//	  const int p = cc[icc].p; // degree
	//	  const int n = cc[icc].n; // number of control points
	//	  const double *t = cc[icc].t; // knot values
	//	  yvec[ix] += freeterm[ix] * closedcurve_potential(p, n, t, &(xvec[xoff[icc]]), tc[ix]); // add freeterm
	//	}

      } else {
	
	l2l(c, nterm, icell); // to children
	
      }
      
      for (int ineigh = 0; ineigh < c[icell].nneigh; ineigh ++) {
	const int jcell = c[icell].neighbors[ineigh];
	const int jflag = imleaf(c[jcell]);
	if (iflag == QUADTREE_TRUE || jflag == QUADTREE_TRUE) {
	  for (int kx = c[icell].chead; kx <= c[icell].ctail; kx ++) {
	    const int ix = carray[kx]; // collocation point in icell
	    //	    const int icc = ccc[ix]; // closec curve of collocation point
	    const int icc = - 1; // dummy index to avoid the (near-)singular integral
	    for (int kb = c[jcell].bhead; kb <= c[jcell].btail; kb ++) {
	      const int j = bj[barray[kb]]; // basis in jcell
	      const int jcc = bcc[barray[kb]]; // closed curve of basis
	      const int p = cc[jcc].p; // degree
	      const int n = cc[jcc].n; // number of control points
	      const double *t = cc[jcc].t; // knot values
	      const double2 *cp = cc[jcc].cp; // control points
	      double dval;
	      //	      isolap2d_direct(eps, g, p, n, t, cp, tc, xc, icc, jcc, ix, j, &dval);
	      isolap2d_direct(eps, g, p, n, t, cp, tc, xin, icc, jcc, ix, j, &dval); // tc is never used
	      yvec[ix] += dval * xvec[xoff[jcc] + j % (n - p)]; // keep the sign of dval
	    }
	  }
	}
      }

      
    } // icell
  } // level

}


static void matvec_fmm_incal(const int nin, const int nterm, const int maxlev, const int minlev, const int ncell,
			     const int *barray, const int *carray, const int *levsta, const int *levend, Cell *c, 
			     const double eps, const GaussOne *g,
			     const Closedcurve *cc, const int *bcc, const int *bj, const int *xoff,
			     const double *tc, const double2 *xin, const double *xvec, double *yvec)
{
  for (int icell = 0; icell < ncell; icell ++) {
    c[icell].mcoe = (dcomplex *)calloc(nterm, sizeof(dcomplex)); // initialise
    c[icell].lcoe = (dcomplex *)calloc(nterm, sizeof(dcomplex)); // initialise
  }
  
  //  for(int ix = 0; ix < ncol; ix ++) {
  for(int ix = 0; ix < nin; ix ++) {
    yvec[ix] = 0; // initialise
  }
  
  ////////////////////////////////////
  MESG("Enter upward\n");
  ////////////////////////////////////
  upward(nterm, maxlev, minlev,
	 barray, cc, bcc, bj, xoff,
	 levsta, levend, c,
	 g, xvec);
  ////////////////////////////////////
  //  MESG("Enter downward\n");
  ////////////////////////////////////
  //  downward(nterm, maxlev, minlev,
  //	   barray, carray, levsta, levend, c, 
  //	   eps, g,
  //	   cc, bcc, bj, xoff, ccc,
  //	   tc, xc, freeterm, xvec, yvec);
  ////////////////////////////////////
  MESG("Enter downward_incal\n");
  ////////////////////////////////////
  downward_incal(nterm, maxlev, minlev,
		 barray, carray, levsta, levend, c, 
		 eps, g,
		 cc, bcc, bj, xoff,
		 tc, xin, xvec, yvec);

  for(int ix = 0; ix < nin; ix ++) {
    yvec[ix] = - yvec[ix]; // change the sign before adding the external field
  }


  for (int icell = 0; icell < ncell; icell ++) {
    free(c[icell].mcoe);
    free(c[icell].lcoe);
  }
}


void isolap2d_incal_fmm(Isolap2d *p, const int nin, const double2 *xin, const double *xvec, double *yvec)
{
  Params *pa = p->params;
  Tree *tree = p->tree;

  matvec_fmm_incal(nin, pa->nterm, tree->maxlev, tree->minlev, tree->ncell, // N(=ncol) => nin
		   tree->barray, tree->carray, tree->levsta, tree->levend, tree->cell, 
		   pa->eps, p->gauss, 
		   p->cc, p->bcc, p->bj, p->xoff,
		   p->tc, xin, xvec, yvec); // p->xc => xin, p->freeterm => removed

}
