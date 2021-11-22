#include "isolap2d.h"


int isolap2d_comp_num_collocation_points(const int ncc, const Closedcurve *cc)
{
  int ncol = 0; // Initialise the number of collocation points (or unknowns)

  for (int icc = 0; icc < ncc; icc ++) {
    const int p = cc[icc].p;
    const int n = cc[icc].n;
    ncol += n - p;
  }
  return ncol;
}


void isolap2d_comp_collocation_points(const int ncc, const Closedcurve *cc, double *tc, double2 *xc, int *ccc, int *xoff)
{
  /* Create n-p collocation points for every closed curve */
  
  int ncol = 0;
  
  /* Loop over closed curves */
  for (int icc = 0; icc < ncc; icc ++) {
    
    const int p = cc[icc].p;
    const int n = cc[icc].n;
    const double *t = cc[icc].t;
    const double2 *cp = cc[icc].cp;

    xoff[icc] = ncol; // offset to xvec

    /* Loop over collocation points */
    for (int i = 0; i < n - p; i ++) {
      
      //    double x = 0.0;
      //    for (int j = 0; j < p; j ++) {
      //      x += t[(i + p) + j]; // skip first p knots
      //    }
      //    x /= p;
      //    tc[i] = x;
      //    xc[i] = bspline_position(p, n, t, cp, x);
      
      tc[ncol] = t[p + i]; // placed at a knot
      xc[ncol] = closedcurve_position(p, n, t, cp, tc[ncol]);
      
      ccc[ncol] = icc;

      ncol ++;
      
    }
  }

  ///////////////////////////////////////////////////////////
  ASSERT(ncol == isolap2d_comp_num_collocation_points(ncc, cc));
  ///////////////////////////////////////////////////////////
}


void isolap2d_check_collocation_points(const int ncol, const double2 *xc)
{
  FILE *fp = fopen("collocation.txt", "w");
  FILE *fp2 = fopen("collocation.label", "w");

  for (int i = 0; i < ncol; i ++) {
    fprintf(fp, "%f %f 0.0\n", xc[i].x, xc[i].y);
    fprintf(fp2, "set label '%d' at %f, %f, %f \n", i, xc[i].x, xc[i].y, 0.0 * i);
  }

  fclose(fp);
  fclose(fp2);
  MSG("Created collocation.txt and collocation.label\n");
}


int isolap2d_comp_num_basis_center(const int ncc, const Closedcurve *cc)
{
  int nbasis = 0; // Initialise the number of bases
  for (int icc = 0; icc < ncc; icc ++) {
    nbasis += cc[icc].n;
  }
  return nbasis;
}


//140318void isolap2d_comp_basis_center(const int ncc, const Closedcurve *cc, double2 *xb, int *bcc, int *bj)
void isolap2d_comp_basis_center(Isolap2d *p)
{
  const int ncc = p->ncc;
  const Closedcurve *cc = p->cc;

  const int nbasis = isolap2d_comp_num_basis_center(ncc, cc);
  INFO("nbasis=%d\n", nbasis);

  p->xb = (double2 *)malloc(nbasis * sizeof(double2)); // xb[0:nbasis-1]
  p->bcc = (int *)malloc(nbasis * sizeof(int)); // bcc[0:nbasis-1]
  p->bj = (int *)malloc(nbasis * sizeof(int)); // bj[0:nbasis-1]

  double2 *xb = p->xb;
  int *bcc = p->bcc;
  int *bj = p->bj;

  //140318  int nbasis = 0;
  int ib = 0; // global index of bases

  /* Loop over closed curves */
  for (int icc = 0; icc < ncc; icc ++) {
    
    const int p = cc[icc].p;
    const int n = cc[icc].n;
    const double *t = cc[icc].t;
    const double2 *cp = cc[icc].cp;
    
    /* Loop over bases */
    for (int j = 0; j < n; j ++) {
      
      /* Obtain the parameter that corresponds to the centre of this
	 basis */
      const double tnow = (t[j] + t[j + p + 1]) / 2.0;

      /* Compute the centre */
      //140318      xb[nbasis] = closedcurve_position(p, n, t, cp, tnow);
      xb[ib] = closedcurve_position(p, n, t, cp, tnow);

      /* Store icc and j */
      //140318      bcc[nbasis] = icc;
      //140318      bj[nbasis] = j;
      bcc[ib] = icc;
      bj[ib] = j;

      /* Update */
      //140318      nbasis ++;
      ib ++;

    }
  }

  ASSERT(ib == nbasis);

}


void isolap2d_free_basis_center(Isolap2d *p)
{
  free(p->xb);
  free(p->bcc);
  free(p->bj);
}


//140318void isolap2d_check_basis_center(const int nbasis, const double2 *xb)
void isolap2d_check_basis_center(const Isolap2d *p)
{
  FILE *fp = fopen("basis.txt", "w");
  FILE *fp2 = fopen("basis.label", "w");

  const int nbasis = isolap2d_comp_num_basis_center(p->ncc, p->cc);
  const double2 *xb = p->xb;

  for (int i = 0; i < nbasis; i ++) {
    fprintf(fp, "%f %f 0.0\n", xb[i].x, xb[i].y);
    fprintf(fp2, "set label 'B%d' at %f, %f, %f \n", i, xb[i].x, xb[i].y, 0.0);
  }
  fclose(fp);
  fclose(fp2);

  INFO("Created basis.txt and basis.label\n");
}


void isolap2d_direct(const double eps, const GaussOne *g, const int p, const int n, const double *t, const double2 *cp, const double *tc, const double2 *xc, const int icc, const int jcc, const int ix, const int j, double *dval)
{
  /* Compute the coeffecient of u_j (j=0...n-1) for a given
     collocation point x */

  /* Constant */
  const double eps2 = eps * eps;

  /* Compute the collocation point */
  const double2 x = xc[ix];

  /* Set up for Guassian quadrature */
  const int ng = g->n;
  const double *gx = g->x;
  const double *gw = g->w;

  /* Initialise */
  *dval = 0.0;

  /* Loop over the sub-intervals in the support of b_{j,p}, that is,
     [t_j,t_{j+p+1}) */
  for (int k = 0; k <= p; k ++) {

    /* Choose the sub-interval [t_{j+k},t_{j+k+1}) in [t_p,t_n) */
    if (j + k >= p && j + k + 1 <= n) {
      
      /* Jacobian wrt Gaussian quadrature (this is not Jacobian wrt
	 the curve parameter) */
      const double Jg = (t[j + k + 1] - t[j + k]) / 2.0;
      
      /* Loop over Gaussian abscissa */
      for (int ig = 0; ig < ng; ig ++) {
	
	/* Parameter corresponding to the ig th abscissa, which is in
	   [t_{j+k},t_{j+k+1}) */
	const double tnow = ((1.0 - gx[ig]) * t[j + k] + (1.0 + gx[ig]) * t[j + k + 1]) / 2.0;
	
	/* Integral point y and its derivative dy/dt */
	const double2 y = closedcurve_position(p, n, t, cp, tnow);
	const double2 dy = closedcurve_position_derivative(p, n, t, cp, tnow);
	
	/* Compute r:=|x-y| */
	const double2 vr = sub2(x, y);
	const double rr = dot2(vr, vr);
	
	/* Compute basis b_{j,p} */
	const double bj = closedcurve_basis(p, n, t, j, tnow);

	/* Compute kernel */
	double kern;
	if (rr < eps2 && icc == jcc) { // if y is close to x and they are in the same closed curve

	  const double2 ddy = closedcurve_position_derivative2(p, n, t, cp, tnow);
	  const double2 dddy = closedcurve_position_derivative3(p, n, t, cp, tnow);
	  kern = ((- dy.x * ddy.y + dy.y * ddy.x) / 2.0 + (2.0 * (- dy.x * dddy.y + dy.y * dddy.x)) / 6.0 * (tnow - tc[ix])) / ((dy.x * dy.x + dy.y * dy.y) + (dy.x * ddy.x + dy.y * ddy.y) * (tnow - tc[ix])) * bj;
	} else {

	  kern = (vr.x * dy.y - vr.y * dy.x) / rr * bj;

	}
	
	/* Accumulate */
	*dval += kern * gw[ig] * Jg;
	
      } // ig

    }

  } // k

  *dval /= 2.0 * PI;

}


static void matvec(const double eps, const GaussOne *g, const int ncc, const Closedcurve *cc, const int ncol, const double *tc, const double2 *xc, const double *freeterm, const double *xvec, double *yvec)
{
  /* Execute the matrix-vector product in terms of Neumann problem,
     where the normal derivative of potential is assumed to be zero on
     the boundary */
  
  /* Initialise */
  for(int ix = 0; ix < ncol; ix ++) {
    yvec[ix] = 0.0;
  }
  
  int ioff = 0; // offset for collocation points

  /* Loop over closed curves */
  for (int icc = 0; icc < ncc; icc ++) {

    const int pi = cc[icc].p; // degree
    const int ni = cc[icc].n; // number of control points
    const double *ti = cc[icc].t; // knot values

    /* Loop over collocation points */
#if defined(_OPENMP)
#pragma omp parallel for // OpenMP DEFINED LOOP WAS PARALLELIZED.
#endif
    for (int i = 0; i < ni - pi; i ++) {

      const int ix = ioff + i; // 0<=ix<ncol
      
      int joff = 0; // offset for basis functions
    
      /* Loop over closed curves */
      for (int jcc = 0; jcc < ncc; jcc ++) {
	
	const int p = cc[jcc].p; // degree
	const int n = cc[jcc].n; // number of control points
	const double *t = cc[jcc].t; // knot values
	const double2 *cp = cc[jcc].cp; // control points
	
	/* Loop over base functions */
	for (int j = 0; j < n; j ++) {
	  
	  /* Compute the interaction */
	  double dval;
	  isolap2d_direct(eps, g, p, n, t, cp, tc, xc, icc, jcc, ix, j, &dval);
	  
	  /* Accumulation */
	  yvec[ix] += dval * xvec[joff + j % (n - p)]; // j%(n-p) in [0,n-p)
	  
	} // j
	
	joff += n - p;
	
      } // jcc
      
      /* Add free term */
      yvec[ix] += freeterm[ix] * closedcurve_potential(pi, ni, ti, &(xvec[ioff]), tc[ix]);

    } // i

    ioff += ni - pi;

  } // icc

}


void isolap2d_matvec(const int N, const void *A, const double *xvec, double *yvec)
{
  Isolap2d *p = (Isolap2d *)A;

#if defined(ORIGINAL)
  matvec(p->params->eps, p->gauss, p->ncc, p->cc, N, p->tc, p->xc, p->freeterm, xvec, yvec);
#else
  matvec_nbie(p->gauss, p->ncc, p->cc, N, p->tc, p->xc, p->freeterm, xvec, yvec);
#endif

}


void isolap2d_comp_freeterm(Isolap2d *p, const int ncol)
{
#if !defined(SET_FREE_TERM)

  /* Compute the free term thorough matvec */

  double *mone = (double *)malloc(ncol * sizeof(double)); // set minus one as xvec
  double *ytmp = (double *)malloc(ncol * sizeof(double)); // store free term temporalily; initilised by matvec

  for (int ix = 0; ix < ncol; ix ++) {
    mone[ix] = - 1.0;
    p->freeterm[ix] = 0.0; // initialise
  }

#if defined(FMM)
  isolap2d_matvec_fmm(ncol, (const void *)p, mone, ytmp);
#else
  isolap2d_matvec(ncol, (const void *)p, mone, ytmp);
#endif

  if (p->params->ninc != NULL_UINC) {
    for (int ix = 0; ix < ncol; ix ++) {
      p->freeterm[ix] = ytmp[ix] + 1.0; // contribution from external field
    }
  }

  free(mone);
  free(ytmp);

#else // set 1/2

  for (int ix = 0; ix < ncol; ix ++) {
    p->freeterm[ix] = 0.5; // smooth
  }

#endif

  ///////////////////////////////////////////////////////////
#if defined(MYDEBUG)
  DBG("free term ncol=%d:\n", ncol);
  for (int ix = 0; ix < ncol; ix ++) {
    //    DBG("ix=%d: %f\n", ix, p->freeterm[ix]);
    DBG("ix=%d: %lf %e\n", ix, p->freeterm[ix], p->freeterm[ix] - 0.5); // difference
  }
#endif
  ///////////////////////////////////////////////////////////
}



int isolap2d_comp_num_element_center(const int ncc, const Closedcurve *cc)
{
  int nelem = 0; // Initialise the number of isogeometric elements
  for (int icc = 0; icc < ncc; icc ++) {
    nelem += cc[icc].n - cc[icc].p; // I_0:=[t_p,t_{p+1}), ..., I_{n-p-1}:=[t_{n-1},t_n), where t_n=t_p
  }
  return nelem;
}


void isolap2d_comp_element_center(Isolap2d *p)
{
  const int ncc = p->ncc;
  const Closedcurve *cc = p->cc;

  const int nelem = isolap2d_comp_num_element_center(ncc, cc);
  INFO("nelem=%d\n", nelem);

  p->xb = (double2 *)malloc(nelem * sizeof(double2)); // xb[0:nelem-1]
  p->bcc = (int *)malloc(nelem * sizeof(int)); // bcc[0:nelem-1]
  p->bj = (int *)malloc(nelem * sizeof(int)); // bj[0:nelem-1]

  double2 *xb = p->xb;
  int *bcc = p->bcc;
  int *bj = p->bj;

  int ie = 0; // global index of elements

  /* Loop over closed curves */
  for (int icc = 0; icc < ncc; icc ++) {
    
    const int p = cc[icc].p;
    const int n = cc[icc].n;
    const double *t = cc[icc].t;
    const double2 *cp = cc[icc].cp;
    
    /* Loop over elements (sections) in this curve */
    for (int j = 0; j < n - p; j ++) {
      
      /* Obtain the parameter that corresponds to the centre of this element */
      const double tnow = (t[j + p] + t[j + p + 1]) / 2.0;

      /* Compute the centre */
      xb[ie] = closedcurve_position(p, n, t, cp, tnow);

      /* Store icc and j */
      bcc[ie] = icc;
      bj[ie] = j;

      /* Update */
      ie ++;

    }
  }

  ASSERT(ie == nelem);

}


void isolap2d_free_element_center(Isolap2d *p)
{
  free(p->xb);
  free(p->bcc);
  free(p->bj);
}


void isolap2d_check_element_center(const Isolap2d *p)
{
  FILE *fp = fopen("elem.txt", "w");
  FILE *fp2 = fopen("elem.label", "w");

  const int nelem = isolap2d_comp_num_element_center(p->ncc, p->cc);
  const double2 *xb = p->xb;

  for (int i = 0; i < nelem; i ++) {
    fprintf(fp, "%f %f 0.0\n", xb[i].x, xb[i].y);
    fprintf(fp2, "set label 'B%d' at %f, %f, %f \n", i, xb[i].x, xb[i].y, 0.0);
  }
  fclose(fp);
  fclose(fp2);

  INFO("Created elem.txt and elem.label\n");
}
