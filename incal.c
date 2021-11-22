#include "isolap2d.h"

static void incal(const int ninc, const double eps, const GaussOne *g, const int ncc, const Closedcurve *cc,
		  const double *tc, const double2 *xc, const int nin, const double2 *xin, const double *xvec, double *yvec)
{
  /*
    Computes the potentials at internal points
  */

  /* Initialise */
  for (int i = 0; i < nin; i ++) {
    yvec[i] = 0.0;
  }
  
  /* Set dummy index for observational closed curve */
  const int icc = - 1;

  /* Loop over internal points */
#if defined(_OPENMP)
#pragma omp parallel for // OpenMP DEFINED LOOP WAS PARALLELIZED.
#endif
  for (int ix = 0; ix < nin; ix ++) {
    
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
	isolap2d_direct(eps, g, p, n, t, cp, tc, xin, icc, jcc, ix, j, &dval); // xc is replaced with xin; tc is never used
	  
	/* Accumulation (minus stands for moving DLP to RHS) */
	yvec[ix] -= dval * xvec[joff + j % (n - p)]; // j%(n-p) in [0,n-p)
	  
      } // j
	
      joff += n - p;
	
    } // jcc

  } // ix

}


void isolap2d_incal(const Isolap2d *p, const int nin, const double2 *xin, const double *xvec, double *yvec)
{
  Params *pa = p->params;

  incal(pa->ninc, pa->eps, p->gauss, p->ncc, p->cc, p->tc, p->xc, nin, xin, xvec, yvec);
}


void isolap2d_incal_input(const Isolap2d *p, const char *infile, int *nin, double2 **xin)
{
  FILE *fp = fopen(infile, "r");
  fscanf(fp, "%d", nin);
  INFO("nin=%d\n", *nin);

  *xin = (double2 *)malloc(*nin * sizeof(double2));

  for (int i = 0; i < *nin; i ++) {
    int num;
    double tmpx, tmpy;
#if defined(DISABLE_INCAL_READ_BOUNDARY_MARKER)
    fscanf(fp, "%d %lf %lf", &num, &tmpx, &tmpy);
#else
    int id; // boundary marker
    fscanf(fp, "%d %lf %lf %d", &num, &tmpx, &tmpy, &id);
#endif
    (*xin)[i].x = tmpx;
    (*xin)[i].y = tmpy;
  }

  fclose(fp);
}


void isolap2d_incal_input_check(const int nin, const double2 *xin)
{
  FILE *fp = fopen("internalpoints.txt", "w");
  FILE *fp2 = fopen("internalpoints.label", "w");
  for (int i = 0; i < nin; i ++) {
    fprintf(fp, "%f %f 0.0\n", xin[i].x, xin[i].y);
    fprintf(fp2, "set label 'I%d' at %f, %f, 0.0\n", i, xin[i].x, xin[i].y);
  }
  fclose(fp);
  fclose(fp2);
  MSG("Created internalpoints.txt and internalpoints.label\n");
}


void isolap2d_incal_output(const int nin, const double2 *xin, const double *yvec)
{
  for (int ix = 0; ix < nin; ix ++) {
    printf("%d %24.15e %24.15e %15.7e %15.7e #INCAL\n", ix, yvec[ix], 0.0, xin[ix].x, xin[ix].y); // qval is zero due to B.C.
  }
}
