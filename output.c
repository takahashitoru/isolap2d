#include "isolap2d.h"

void isolap2d_output(const int ncc, const Closedcurve *cc, const double *tc, const double2 *xc, const double *xvec)
{
  int ioff = 0; // offset for collocation points
  int ix = 0;   // counter for collocation points  

  /* Loop over closed curve */
  for (int icc = 0; icc < ncc; icc ++) {
    
    const int p = cc[icc].p;
    const int n = cc[icc].n;
    const double *t = cc[icc].t;
  
    /* Loop over collocation points */
    for (int i = 0; i < n - p; i ++) {

      /* Potential at the current collocation point */
      const double uval = closedcurve_potential(p, n, t, &(xvec[ioff]), tc[ix]);

      /* Print the potential */
      printf("%d %24.15e %24.15e %15.7e %15.7e\n", ix, uval, 0.0, xc[ix].x, xc[ix].y); // qval is zero due to B.C.            

      ix ++;

    } // i

    ioff += n - p;

  } // icc

}


void isolap2d_error(const int ncc, const Closedcurve *cc, const double *tc, const double2 *xc, const double *xvec)
{
  /* Compute the relative L2-error and the maximum error for the
     external problem where 'x' is supplied as the uniform external
     field (strength is one) to the circular boundary with any radius
     at the origin */

  double err = 0.0;
  double sum = 0.0;
  double errmax = 0.0;
  double uanamax = 0.0;
  
  int ioff = 0; // offset for collocation points
  int ix = 0;   // counter for collocation points

  /* Loop over closed curve */
  for (int icc = 0; icc < ncc; icc ++) {
    
    const int p = cc[icc].p;
    const int n = cc[icc].n;
    const double *t = cc[icc].t;

    /* Loop over collocation points */
    for (int i = 0; i < n - p; i ++) {
    
      /* Potential at the current collocation point */
      const double uval = closedcurve_potential(p, n, t, &(xvec[ioff]), tc[ix]);
      
      /* Analytical solution */
      //const double uana = 2.0 * xc[ix].x;
      const double uana = 2.0 * xc[ix].y;

    
      /* Stack */
      err += (uval - uana) * (uval - uana);
      sum += uana * uana;
      
      /* Maximum norm */
      errmax = MAX(errmax, (uval - uana) * (uval- uana));
      uanamax = MAX(uanamax, uana * uana);
    
      ix ++;

    } // i

    ioff += n - p;

  } // icc
  
  /* Result */
  INFO("L2error=%24.15e\n", sqrt(err / sum));
  INFO("Maxerror=%24.15e\n", sqrt(errmax / uanamax));
  
}
