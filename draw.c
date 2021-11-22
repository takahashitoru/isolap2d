#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "isolap2d.h"

// gcc -std=c99 -I. draw.c closedcurve.c double2.c -lm

int main()
{
  int p = 2; // order of B-spline basis functions (quadratic if 2)

#if(1)
  int n = 17 + p; // 'C' shape
#endif
#if(0)
  int n = 9; // number of control points (0 to n-1)
#endif
#if(0)
  int n = 10000;
#endif

  fprintf(stderr, "p=%d\n", p);
  fprintf(stderr, "n=%d\n", n);
  
  //  /* Compute the knot vector */
  //
  //  double *t = bspline_knot_alloc(p, n);
  //  bspline_knot_comp(p, n, t);

  /* Define n control points; place the last p points on the first p
     points */
  
  double2 *cp = (double2 *)malloc(n * sizeof(double2));

  // first n-p points
#if(1)
  cp[ 0].x = -191; cp[ 0].y =  -76;
  cp[ 1].x =  -77; cp[ 1].y = -198;
  cp[ 2].x =   92; cp[ 2].y = -198;
  cp[ 3].x =  207; cp[ 3].y =  -36;
  cp[ 4].x =  123; cp[ 4].y =  -11;
  cp[ 5].x =  112; cp[ 5].y =  -73;
  cp[ 6].x =   25; cp[ 6].y = -146;
  cp[ 7].x =  -99; cp[ 7].y = -112;

  cp[ 8].x = -144; cp[ 8].y =    0;

  cp[ 9].x =  -99; cp[ 9].y =  112;
  cp[10].x =   25; cp[10].y =  146;
  cp[11].x =  112; cp[11].y =   73;
  cp[12].x =  123; cp[12].y =   11;
  cp[13].x =  207; cp[13].y =   36;
  cp[14].x =   92; cp[14].y =  198;
  cp[15].x =  -77; cp[15].y =  198;
  cp[16].x = -191; cp[16].y =   96;

  //	248,448,134,346,126,174,273,45,455,96,523,231,448,248,437,185,
  //	350,117,229,145,181,255,226,362,350,396,437,323,448,261,532,286,
  //	417,448,248,448
#endif
#if(0)
  cp[0].x = 0; cp[0].y = 0;
  cp[1].x = 0; cp[1].y = 2;
  cp[2].x = 2; cp[2].y = 3;
  cp[3].x = 4; cp[3].y = 0;
  cp[4].x = 6; cp[4].y = 3;
  cp[5].x = 8; cp[5].y = 2;
  cp[6].x = 8; cp[6].y = 0;
#endif
#if(0)
  for (int i = 0; i < n - p; i ++) {
    double t = 2.0 * PI * i / (n - p - 1);
    //    double a = 1.0 + 0.8 * cos(20 * t);
    double a = 1.0 + 0.4 * RAND(0.0, 1.0) * cos(20 * t);
    cp[i].x = a * cos(t);
    cp[i].y = a * sin(t);
  }
#endif
#if(0)
  for (int i = 0; i < n - p; i ++) {
    cp[i].x = RAND(- 1.0, 1.0);
    cp[i].y = RAND(- 1.0, 1.0);
  }
#endif

  for (int i = n - p; i < n; i ++) { // cp[n-p:n-1]
    cp[i] = cp[i - (n - p)]; // wrap the last p points on the first p points
  }

  /* Set a closed curve */
  Closedcurve *cc = closedcurve_allocate(1);
  cc->p = p;
  cc->n = n;
  cc->cp = closedcurve_cp_alloc(n);  
  for (int i = 0; i < n; i ++) {
    (cc->cp)[i].x = cp[i].x;
    (cc->cp)[i].y = cp[i].y;
  }
  cc->t = closedcurve_knot_alloc(p, n);
  closedcurve_knot_comp(p, n, cc->t);



  /* Check control points */
  FILE *fp = fopen("cp.txt", "w");
  FILE *fp2 = fopen("cp.label", "w");
  for (int i = 0; i < n; i ++) { // cp[0:n-1]
    fprintf(stderr, "i=%d cp=%f %f 0.0\n", i, cp[i].x, cp[i].y);
    fprintf(fp, "%f %f 0.0\n", cp[i].x, cp[i].y);
    fprintf(fp2, "set label 'C%d' at %f, %f, %f \n", i, cp[i].x, cp[i].y, 0.1 * i);
  }
  fclose(fp);
  fclose(fp2);


  /* Check knots */

  closedcurve_knot_check(1, cc);

  /* Draw the boundary */

  closedcurve_draw_boundary(1, cc);
  
  closedcurve_cp_free(cc->cp);
  closedcurve_knot_free(cc->t);
  closedcurve_free(cc);

  free(cp);

  return(EXIT_SUCCESS);
}
