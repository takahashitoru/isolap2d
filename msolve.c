#include "isolap2d.h"

//LATER static void msolve_point_jacobi(const int nb, const double eps, const int *nod, const int *ibc, const double2 *x, const double2 *xg, const double *rvec, double *zvec)
//LATER {
//LATER   for (int ix = 0; ix < nb; ix ++) {
//LATER 
//LATER     double2 xloc[2], vn;
//LATER     double length;
//LATER     ldedat(nod, x, ix, xloc, &vn, &length); // iy=ix
//LATER     
//LATER     double sval, dval;
//LATER     direct(eps, xloc, vn, xg[ix], &sval, &dval);
//LATER 
//LATER     double uval, qval;
//LATER     bval(ibc, xg, ix, vn, BVAL_PRECOND, NULL, &uval, &qval);
//LATER 
//LATER     zvec[ix] = rvec[ix] / (dval * uval - sval * qval);
//LATER 
//LATER   }
//LATER }

#if(0)
static void msolve_point_jacobi()
{
  /* This is not evident. Which collocation point relates to a certain boundary data u_j? */
}
#endif

/* Solve Mz=r for z */

void isolap2d_msolve(const int N, const void *M, const double *r, double *z)
{
  Isolap2d *p = (Isolap2d *)M;
  Params *pa = p->params;

  if (pa->jpre == 0) {

    for (int ix = 0; ix < N; ix++) { // do nothing, that is, let M=I
      z[ix] = r[ix];
    }

    //LATER   } else if (pa->jpre == 1) {
    //LATER     msolve_point_jacobi(p->nb, pa->eps, p->nod, p->ibc, p->x, p->xg, r, z);

  } else {

    MESG("Unregistered jpre was used. Exit.\n");
    exit(EXIT_FAILURE);

  }
}
