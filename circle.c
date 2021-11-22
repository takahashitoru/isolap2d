#include "isolap2d.h"

static void create_circle(const int p, const int n, const double radius, const double cx, const double cy, double2 *cp)
{
  /* Place n control points on a circle. Note that the resulting
     clolsed curve will have a smaller radius than the specified one
     because B-spline functions are not capable of reprsenting an
     exact circle */

  for (int i = 0; i < n - p; i ++) { // cp[0:n-p-1]
    double theta = - 2.0 * PI * i / (n - p); // minus for external problem
    cp[i].x = radius * cos(theta) + cx;
    cp[i].y = radius * sin(theta) + cy;
  }  

  for (int i = n - p; i < n; i ++) { // cp[n-p:n-1]
    cp[i] = cp[i - (n - p)]; // wrap the last p points on the first p points
  }

}

//120124void isolap2d_input(Isolap2d *p)
void isolap2d_input_circle(Isolap2d *p)
{
  /* Define the number of closed curves in the current geometry */
  p->ncc = CIRCLE_NUM_CLOSED_CURVES;
  INFO("ncc=%d\n", p->ncc);

  /* Allocate closed curves */
  p->cc = closedcurve_allocate(p->ncc);

  /* Loop over closed curves */
  for (int icc = 0; icc < p->ncc; icc ++) {

    /* Define the degree of B-spline functions */
    (p->cc)[icc].p = CIRCLE_SPLINE_DEGREE;

    /* Define the number of control points */
    (p->cc)[icc].n = CIRCLE_NUM_CONTROL_POINTS;

    /* Print infomation */
    INFO("icc=%d: p=%d n=%d\n", icc, (p->cc)[icc].p, (p->cc)[icc].n);

    /* Allocate n control points */
    (p->cc)[icc].cp = closedcurve_cp_alloc((p->cc)[icc].n);

    /* Define the radius and center */
    const double radius = CIRCLE_RADIUS;
    const double cx = icc % 10;
    const double cy = icc / 10;

    /* Create n control points */
    create_circle((p->cc)[icc].p, (p->cc)[icc].n, radius, cx, cy, (p->cc)[icc].cp);

  }

}

void isolap2d_input_circle_info()
{
  INFO("CIRCLE_NUM_CLOSED_CURVES=%d\n", CIRCLE_NUM_CLOSED_CURVES);
  INFO("CIRCLE_RADIUS=%f\n", (double)CIRCLE_RADIUS);
  INFO("CIRCLE_SPLINE_DEGREE=%d\n", CIRCLE_SPLINE_DEGREE);
  INFO("CIRCLE_NUM_CONTROL_POINTS=%d\n", CIRCLE_NUM_CONTROL_POINTS);
}
