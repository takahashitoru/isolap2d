//
// dgmres.c
//
// This is a preconditioned GMRES(m) code in double precision.
// This was obtained from gmres.h in cpptemplates.
//
// 2009. 2.25: nmv and nms were added.
// 2008. 6. 3: The definition of nullx0 was changed.
// 2008. 2.29: A new parameter nullx0 was added.
// 2008. 2.26: A bug was fixed. msolve is available.
// 2008. 2.25: static was added to the internal functions.
// 2008. 2.23: obtained from gmres.h in cpptemplates.
//  
// Note: Following is the preface of gmres.h
//*****************************************************************
// Iterative template routine -- GMRES
//
// GMRES solves the unsymmetric linear system Ax = b using the 
// Generalized Minimum Residual method
//
// GMRES follows the algorithm described on p. 20 of the 
// SIAM Templates book.
//
// The return value indicates convergence within max_iter (input)
// iterations (0), or no convergence within max_iter iterations (1).
//
// Upon successful return, output arguments have the following values:
//  
//        x  --  approximate solution to Ax = b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//  
//*****************************************************************

#include "dgmres.h"


void dgmres_info()
{
#if defined(DGMRES_RATE_CONVERGENCE)
  fprintf(stderr, "# %s: DGMRES_RATE_CONVERGENCE\n", __FILE__);
  //140225#endif
#if defined(DGMRES_KEEP_RESIDUALS)
  fprintf(stderr, "# %s: DGMRES_KEEP_RESIDUALS = %d\n", __FILE__, DGMRES_KEEP_RESIDUALS);
#endif
#endif
}


static double dot(const int n, const double *x, const double *y)
{
  double s = 0.0;
  for (int i = 0; i < n; i ++) {
    s += x[i] * y[i];
  }
  return s;
}


static double norm(const int n, const double *x)
{
  return sqrt(dot(n, x, x));
}


static void Update(double *x, const int k, const double *H, const double *s, const double *v, const int m, const int n)
{
  double *y = (double *)malloc((k + 1) * sizeof(double));
  for (int i = 0; i <= k; i ++) {
    y[i] = s[i];
  }

  /* Backsolve */
  for (int i = k; i >= 0; i --) {
    y[i] /= H[i * m + i];
    for (int j = i - 1; j >= 0; j --) {
      y[j] -= H[j * m + i] * y[i];
    }
  }

  for (int j = 0; j <= k; j ++) {
    for (int l = 0; l < n; l ++) {
      x[l] += v[j * n + l] * y[j];
    }
  }
  free(y);
}


static void ApplyPlaneRotation(double *dx, double *dy, const double cs, const double sn)
{
  double temp = cs * (*dx) + sn * (*dy);
  *dy = -sn * (*dx) + cs * (*dy);
  *dx = temp;
}


static void GeneratePlaneRotation(const double dx, const double dy, double *cs, double *sn)
{
  if (dy == 0.0) {
    *cs = 1.0;
    *sn = 0.0;
  } else if (fabs(dy) > fabs(dx)) {
    double temp = dx / dy;
    *sn = 1.0 / sqrt(1.0 + temp * temp);
    *cs = temp * (*sn);
  } else {
    double temp = dy / dx;
    *cs = 1.0 / sqrt(1.0 + temp * temp);
    *sn = temp * (*cs);
  }
}


static void messageOnReturn(const int niter, const double resid, const int nmv, const int nms)
{
#ifndef GMRES_SPPRESS_MESSAGE
  fprintf(stderr, "# %s: niter = %d\n", __FILE__, niter);
  fprintf(stderr, "# %s: resid = %7.3e\n", __FILE__,  resid);
  fprintf(stderr, "# %s: nmv   = %d\n", __FILE__, nmv);
  fprintf(stderr, "# %s: nms   = %d\n", __FILE__, nms);
#endif
}


static int check_covergence(const double tol, double *resid, const int maxresiduals, double *residuals)
{
#if defined(DGMRES_RATE_CONVERGENCE) // check by the rate of residuals

  for (int i = 0; i < maxresiduals - 1; i ++) { // update
    residuals[i + 1] = residuals[i];
  }
  residuals[0] = *resid;

  if (residuals[maxresiduals - 1] == DGMRES_DUMMY_RESIDUAL) { // the last one is still not updated

    return DGMRES_FALSE; // unconverged

  } else {

    double rate = 0.0;
    for (int i = 0; i < maxresiduals - 1; i ++) {
      rate += (residuals[i + 1] - residuals[i]) / residuals[i + 1]; // residuals[i+1]>residuals[i]
    }
    rate /= (maxresiduals - 1); // averaged rate

    if (rate <= tol) {
      return DGMRES_TRUE; // converged
    } else {
      return DGMRES_FALSE; // unconverged
    }

  }

#else // check by the current residual

  if (*resid <= tol) {
    return DGMRES_TRUE; // converged
  } else {
    return DGMRES_FALSE; // unconverged
  }

#endif

}



int dgmres(const int n, const void *A, const void *M, double *x, const double *b,  
	   const int m, const double tol, const int maxiter, const int nullx0,
	   void (*matvec)(const int, const void *, const double *, double *),  
	   void (*msolve)(const int, const void *, const double *, double *),  
	   double *resid, int *niter, int *nmv, int *nms)
{
  /* Initialise numbers */
  *nmv = 0; // number of performing matvec
  *nms = 0; // number of performing msolve

  /* Allocate memories */
  //120116  double *H = (double *)malloc((m + 1) * m * sizeof(double)); // H[0:m][0:m-1] : Hessenberg matrix
  double *z = (double *)malloc(n * sizeof(double)); // z[0:n-1] : working array
  double *r = (double *)malloc(n * sizeof(double)); // r[0:n-1] : working array (residual vector)

  /*
    Start
  */

  /* Compute the norm of the preconditioned rhs vector */
  msolve(n, M, b, z); (*nms) ++; // z:=M^-1*b (solve z from Mz=b)
  const double normb = norm(n, z); // |M^-1*b|

  if (normb == 0.0) {
    *resid = NAN;
    *niter = 0;
    free(z); free(r);
    return DGMRES_FAIL; // M^-1 or b is wrong
  }

  /* If nullx0 is true, the initial guess x0 is assumed to be zero,
     and then matvec can be skipped */
  if (nullx0 != DGMRES_TRUE) {
    matvec(n, A, x, z); (*nmv) ++; // z:=A*x0
    for (int i = 0; i < n; i ++) {
      z[i] = b[i] - z[i]; // b-A*x0
    }
  } else { // matvec may be skipped if x0=0
    for (int i = 0; i < n; i ++) {
      z[i] = b[i]; // b-A*x0, where x0=0
    }
  }
  msolve(n, M, z, r); (*nms) ++; // Solve r from Mr=b-A*x0
  double beta = norm(n, r); // beta:=|r0|

  //120116  if (normb == 0.0) {
  //120116    normb = 1.0;
  //120116  }

  const int maxresiduals = DGMRES_KEEP_RESIDUALS;
  double residuals[maxresiduals];
  for (int i = 0; i < maxresiduals; i ++) {
    residuals[i] = DGMRES_DUMMY_RESIDUAL; // initialize
  }
  
  *resid = beta / normb;

  //120116  if (*resid <= tol) {
  if (check_covergence(tol, resid, maxresiduals, residuals) == DGMRES_TRUE) {
    *niter = 0;
    messageOnReturn(*niter, *resid, *nmv, *nms);
    free(z); free(r);
    return DGMRES_SUCCESS;
  }

  /* Allocate memories */
  double *H  = (double *)malloc((m + 1) * m * sizeof(double)); // H[0:m][0:m-1] : Hessenberg matrix
  double *v  = (double *)malloc((m + 1) * n * sizeof(double)); // v[0:m][0:n-1]
  double *s  = (double *)malloc((m + 1) * sizeof(double)); // s[0:m]
  double *cs = (double *)malloc(m * sizeof(double)); // cs[0:m-1]
  double *sn = (double *)malloc(m * sizeof(double)); // sn[0:m-1]
  double *w  = (double *)malloc(n * sizeof(double)); // w[0:n-1]

  int j = 1;
  while (j <= maxiter) {
    
    for (int i = 0; i < n; i ++) {
      v[0 * n + i] = r[i] / beta; // v0:=r/beta
    }
    
    for (int i = 0; i < m + 1; i ++) {
      s[i] = 0.0;
    }
    s[0] = beta;
    
    for (int i = 0; i < m && j <= maxiter; i ++, j ++) {

      matvec(n, A, &v[i * n + 0], z); (*nmv) ++;
      msolve(n, M, z, w); (*nms) ++; // w = M.solve(A * v[i]); (solve w from Mw=Av[i])

      for (int k = 0; k <= i; k ++) {
	H[k * m + i] = dot(n, w, &v[k * n + 0]);
	for (int l = 0; l < n; l ++) {
	  w[l] -= H[k * m + i] * v[k * n + l];
	}
      }
      H[(i + 1) * m + i] = norm(n, w);
      for (int l = 0; l < n; l ++) {
	v[(i + 1) * n + l] = w[l] / H[(i + 1) * m + i];
      }

      for (int k = 0; k < i; k ++) {
	ApplyPlaneRotation(&H[k * m + i], &H[(k + 1) * m + i], cs[k], sn[k]);
      }
      
      GeneratePlaneRotation(H[i * m + i], H[(i + 1) * m + i], &cs[i], &sn[i]);
      ApplyPlaneRotation(&H[i * m + i], &H[(i + 1) * m + i], cs[i], sn[i]);
      ApplyPlaneRotation(&s[i], &s[i + 1], cs[i], sn[i]);
      
      *resid = fabs(s[i + 1]) / normb;
      
#ifndef GMRES_SPPRESS_MESSAGE
      fprintf(stderr, "# %s: iter = %i   residual = %7.3e\n", __FILE__, j, *resid);
#endif

      //120116      if (*resid <= tol) {
      if (check_covergence(tol, resid, maxresiduals, residuals) == DGMRES_TRUE) {
        Update(x, i, H, s, v, m, n);
	*niter = j;
        free(H); free(z); free(r); free(v); free(s); free(cs); free(sn); free(w);
	messageOnReturn(*niter, *resid, *nmv, *nms);
        return DGMRES_SUCCESS;
      }

    } // i

    Update(x, m - 1, H, s, v, m, n);
    matvec(n, A, x, z); (*nmv) ++;

    for (int l = 0; l < n; l ++) {
      z[l] = b[l] - z[l];
    }
    msolve(n, M, z, r); (*nms) ++; // solve r from Mr=z(=b-Ax)

    beta = norm(n, r);
    *resid = beta / normb;

    //120116    if (*resid <= tol) {
    if (check_covergence(tol, resid, maxresiduals, residuals) == DGMRES_TRUE) {
      *niter = j;
      free(H); free(z); free(r); free(v); free(s); free(cs); free(sn); free(w);
      messageOnReturn(*niter, *resid, *nmv, *nms);
      return DGMRES_SUCCESS;
    }

  }
  
  free(H); free(z); free(r); free(v); free(s); free(cs); free(sn); free(w);

  messageOnReturn(0, *resid, *nmv, *nms);

  return DGMRES_FAIL;
}
