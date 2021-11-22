#include "isolap2d.h"

static void line(const int ndiv, const double2 center, const double2 vsta, const double2 vend, int *node, double2 *cp)
{
  for (int i = 0; i < ndiv; i ++) {
    double tmp = (double) i / ndiv;
    cp[*node].x = center.x + tmp * (vend.x - vsta.x) + vsta.x;
    cp[*node].y = center.y + tmp * (vend.y - vsta.y) + vsta.y;
    (*node) ++;
  }
}

static void hexa(const int p, const int ndiv, double a, const double2 center, const double scale, int *node, double2 *cp)
{

  double2 v[6];
  double b = a * SQRT3 / 2.0;
  v[0].x = + b;
  v[0].y = + a / 2.0;
  v[1].x = + b;
  v[1].y = - a / 2.0;
  v[2].x = 0.0;
  v[2].y = - b;
  v[3].x = - b;
  v[3].y = - a / 2.0;
  v[4].x = - b;
  v[4].y = + a / 2.0;
  v[5].x = 0.0;
  v[5].y = b;

  for (int i = 0; i < 6; i ++) {
    v[i].x *= scale;
    v[i].y *= scale;
  }

  line(ndiv, center, v[0], v[1], node, cp);
  line(ndiv, center, v[1], v[2], node, cp);
  line(ndiv, center, v[2], v[3], node, cp);
  line(ndiv, center, v[3], v[4], node, cp);
  line(ndiv, center, v[4], v[5], node, cp);
  line(ndiv, center, v[5], v[0], node, cp);

  for (int i = 0; i < p; i ++) { // wrap p points to the first p ones
    double tmp = (double) i / ndiv;
    cp[*node].x = center.x + tmp * (v[1].x - v[0].x) + v[0].x;
    cp[*node].y = center.y + tmp * (v[1].y - v[0].y) + v[0].y;
    (*node) ++;
  }    
  
}

void isolap2d_input_hexa(Isolap2d *p)
{
  const int degree = HEXA_SPLINE_DEGREE;
  const int ndiv = HEXA_NUM_CONTROL_POINTS_PER_EDGE;
  const int num_per_hexa = 6 * ndiv + degree; // number of CP per hexagonal

  const double a = HEXA_RADIUS;
  const double scale = HEXA_SCALE;

  const int nx = HEXA_NUM_FOR_X;
  const int ny = HEXA_NUM_FOR_Y;

  p->ncc = nx * ny; // number of closed curves
  INFO("ncc=%d\n", p->ncc);

  p->cc = closedcurve_allocate(p->ncc);

  int icc = 0; // initialise the index for closed curve

  for (int j = 0; j < ny; j ++) {
    for (int i = 0; i < nx; i ++) {

      (p->cc)[icc].p = degree;
      (p->cc)[icc].n = num_per_hexa;
      //////////////////////////////////////////////////////////////////
      DBG("icc=%d: p=%d n=%d\n", icc, (p->cc)[icc].p, (p->cc)[icc].n);
      //////////////////////////////////////////////////////////////////

      (p->cc)[icc].cp = closedcurve_cp_alloc((p->cc)[icc].n); // control points

      double2 center;
      center.x = SQRT3 * a * i - SQRT3 / 2.0 * a * (j % 2);
      center.y = 1.5 * a * j;

      (p->cc)[icc].x = center; // store the center for the purpose of meshing etc

      int node = 0;

      hexa((p->cc)[icc].p, ndiv, a, center, scale, &node, (p->cc)[icc].cp);
      
      icc ++;
    }
  }
  
}

void isolap2d_input_hexa_info()
{
  INFO("HEXA_SPLINE_DEGREE=%d\n", HEXA_SPLINE_DEGREE);
  INFO("HEXA_NUM_CONTROL_POINTS_PER_EDGE=%d\n", HEXA_NUM_CONTROL_POINTS_PER_EDGE);
  INFO("HEXA_RADIUS=%f\n", HEXA_RADIUS);
  INFO("HEXA_SCALE=%f\n", HEXA_SCALE);
  INFO("HEXA_NUM_FOR_X=%d\n", HEXA_NUM_FOR_X);
  INFO("HEXA_NUM_FOR_Y=%d\n", HEXA_NUM_FOR_Y);
}
