#include "isolap2d.h"

static void cshape(const int p, const double a, const double2 center, const double scale, int *node, double2 *cp)
{
  double2 v[17];
  v[ 0].x = -191; v[ 0].y =  -76;
  v[ 1].x =  -77; v[ 1].y = -198;
  v[ 2].x =   92; v[ 2].y = -198;
  v[ 3].x =  207; v[ 3].y =  -36;
  v[ 4].x =  123; v[ 4].y =  -11;
  v[ 5].x =  112; v[ 5].y =  -73;
  v[ 6].x =   25; v[ 6].y = -146;
  v[ 7].x =  -99; v[ 7].y = -112;

  v[ 8].x = -144; v[ 8].y =    0;

  v[ 9].x =  -99; v[ 9].y =  112;
  v[10].x =   25; v[10].y =  146;
  v[11].x =  112; v[11].y =   73;
  v[12].x =  123; v[12].y =   11;
  v[13].x =  207; v[13].y =   36;
  v[14].x =   92; v[14].y =  198;
  v[15].x =  -77; v[15].y =  198;
  v[16].x = -191; v[16].y =   96;

  const double tmp = scale * a / 400; // default lattice length is 400

  for (int i = 0; i < 17; i ++) { // wrap p points to the first p ones
    cp[*node].x = center.x + tmp * v[i].x;
    cp[*node].y = center.y + tmp * v[i].y;
    (*node) ++;
  }    

  for (int i = 0; i < p; i ++) { // wrap p points to the first p ones
    cp[*node].x = center.x + tmp * v[i].x;
    cp[*node].y = center.y + tmp * v[i].y;
    (*node) ++;
  }    
  
}

void isolap2d_input_cshape(Isolap2d *p)
{
  const int degree = CSHAPE_SPLINE_DEGREE;
  const int num_per_cshape = 17 + degree; // number of CP per c-shaped closed curve

  const double a = CSHAPE_LATTICE;
  const double scale = CSHAPE_SCALE;

  const int nx = CSHAPE_NUM_FOR_X;
  const int ny = CSHAPE_NUM_FOR_Y;

  p->ncc = nx * ny; // number of closed curves
  INFO("ncc=%d\n", p->ncc);

  p->cc = closedcurve_allocate(p->ncc);

  int icc = 0; // initialise the index for closed curve

  for (int j = 0; j < ny; j ++) {
    for (int i = 0; i < nx; i ++) {

      (p->cc)[icc].p = degree;
      (p->cc)[icc].n = num_per_cshape;
      //////////////////////////////////////////////////////////////////
      DBG("icc=%d: p=%d n=%d\n", icc, (p->cc)[icc].p, (p->cc)[icc].n);
      //////////////////////////////////////////////////////////////////

      (p->cc)[icc].cp = closedcurve_cp_alloc((p->cc)[icc].n); // control points

      double2 center;
      //      center.x = SQRT3 * a * i - SQRT3 / 2.0 * a * (j % 2);
      //      center.y = 1.5 * a * j;
      center.x = a * i;
      center.y = a * j;

      (p->cc)[icc].x = center; // store the center for the purpose of meshing etc

      int node = 0;

      cshape((p->cc)[icc].p, a, center, scale, &node, (p->cc)[icc].cp);
      
      icc ++;
    }
  }
  
}

void isolap2d_input_cshape_info()
{
  INFO("CSHAPE_SPLINE_DEGREE=%d\n", CSHAPE_SPLINE_DEGREE);
  INFO("CSHAPE_LATTICE=%f\n", CSHAPE_LATTICE);
  INFO("CSHAPE_SCALE=%f\n", CSHAPE_SCALE);
  INFO("CSHAPE_NUM_FOR_X=%d\n", CSHAPE_NUM_FOR_X);
  INFO("CSHAPE_NUM_FOR_Y=%d\n", CSHAPE_NUM_FOR_Y);
  INFO("CSHAPE_NUM_PER_UNIT=%d\n", CSHAPE_NUM_PER_UNIT);
}
