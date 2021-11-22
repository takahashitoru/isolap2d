#include "isolap2d.h"

void isolap2d_input_check_cp(Isolap2d *p)
{
  /* Check control points */

  FILE *fp = fopen("cp.txt", "w");
  FILE *fp2 = fopen("cp.label", "w");

  for (int icc = 0; icc < p->ncc; icc ++) { 
    for (int i = 0; i < (p->cc)[icc].n; i ++) { // cp[0:n-1]
      fprintf(fp, "%f %f 0.0\n", (p->cc)[icc].cp[i].x, (p->cc)[icc].cp[i].y);
      fprintf(fp2, "set label 'C%d' at %f, %f, %f \n", i, (p->cc)[icc].cp[i].x, (p->cc)[icc].cp[i].y, 0.1 * i);
    }
    fprintf(fp, "\n\n"); // separators
  }

  fclose(fp);
  fclose(fp2);

  MSG("Created cp.txt and cp.label\n");
}


void isolap2d_input_dxf(char *dxffile, Isolap2d *p)
{
  /* Read B-splines from DXF file */
  int nbsp;
  Bspline *bsp = (Bspline *)malloc(DXF_MAX_SPLINES * sizeof(Bspline));
  dxf_read_bspline(dxffile, DXF_MAX_SPLINES, &nbsp, bsp);

  //////////////////////////////////////////////////////////////////////////////////
#if defined(MYDEBUG)
  for (int ibsp = 0; ibsp < nbsp ; ibsp ++) {
    int degree = bsp[ibsp].degree;
    int m = bsp[ibsp].m;
    DBG("ibsp=%d: degree=%d m=%d\n", ibsp, degree, m);
    for (int i = 0; i < m; i ++) {
      DBG("i=%d: x=%f y=%f z=%f\n", i, bsp[ibsp].x[i], bsp[ibsp].y[i], bsp[ibsp].z[i]);
    }
  }
#endif
  //////////////////////////////////////////////////////////////////////////////////

  /* Allocate closed curves */
  p->ncc = nbsp;
  p->cc = closedcurve_allocate(p->ncc); // cc[0:ncc-1]
  Closedcurve *cc = p->cc; // shortcut

  INFO("ncc=%d\n", p->ncc);

  /* Convert data */
  for (int i = 0; i < p->ncc ; i ++) {

    cc[i].p = bsp[i].degree;
    cc[i].n = bsp[i].m + cc[i].p; // increased by p because of wrapping p points

    cc[i].cp = closedcurve_cp_alloc(cc[i].n); // cp[0:n-1]

    for (int j = 0; j < cc[i].n - cc[i].p; j ++) { // cp[0:n-p-1]
      cc[i].cp[j].x = bsp[i].x[j];
      cc[i].cp[j].y = bsp[i].y[j];
    }
    for (int j = cc[i].n - cc[i].p; j < cc[i].n; j ++) { // cp[n-p:n-1]
      cc[i].cp[j] = cc[i].cp[j - (cc[i].n - cc[i].p)]; // wrap the last p points on the first p points
    }

    INFO("icc=%d: p=%d n=%d\n", i, cc[i].p, cc[i].n);
  }

  /* Free */
  dxf_free_bspline(nbsp, bsp);
}
