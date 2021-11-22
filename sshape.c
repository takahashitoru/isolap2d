#include "isolap2d.h"

static void sshape(const int p, const double2 center, int *node, double2 *cp)
{
  for (int i = 0; i < SSHAPE_NDIV; i ++) { // top & inner
    double theta = SSHAPE_TOP_STA_ANGLE + (SSHAPE_TOP_END_ANGLE - SSHAPE_TOP_STA_ANGLE) * i / SSHAPE_NDIV;
    cp[*node].x = center.x + SSHAPE_SCALE * ((1.0 - SSHAPE_WIDTH) * SSHAPE_XRAD * cos(theta));
    cp[*node].y = center.y + SSHAPE_SCALE * ((1.0 - SSHAPE_WIDTH) * SSHAPE_YRAD * sin(theta) + SSHAPE_LATTICE / 4.0);
    (*node) ++;
  }
  for (int i = 0; i < SSHAPE_NDIV; i ++) { // bottom & outer
    double theta = SSHAPE_BOT_STA_ANGLE + (SSHAPE_BOT_END_ANGLE - SSHAPE_BOT_STA_ANGLE) * i / SSHAPE_NDIV;
    cp[*node].x = center.x + SSHAPE_SCALE * ((1.0 + SSHAPE_WIDTH) * SSHAPE_XRAD * cos(theta));
    cp[*node].y = center.y + SSHAPE_SCALE * ((1.0 + SSHAPE_WIDTH) * SSHAPE_YRAD * sin(theta) - SSHAPE_LATTICE / 4.0);
    (*node) ++;
  }
  for (int i = 0; i < SSHAPE_NDIV; i ++) { // bottom & inner
    double theta = SSHAPE_BOT_END_ANGLE + (SSHAPE_BOT_STA_ANGLE - SSHAPE_BOT_END_ANGLE) * i / SSHAPE_NDIV; // reverse
    cp[*node].x = center.x + SSHAPE_SCALE * ((1.0 - SSHAPE_WIDTH) * SSHAPE_XRAD * cos(theta));
    cp[*node].y = center.y + SSHAPE_SCALE * ((1.0 - SSHAPE_WIDTH) * SSHAPE_YRAD * sin(theta) - SSHAPE_LATTICE / 4.0);
    (*node) ++;
  }
  for (int i = 0; i < SSHAPE_NDIV; i ++) { // top & outer
    double theta = SSHAPE_TOP_END_ANGLE + (SSHAPE_TOP_STA_ANGLE - SSHAPE_TOP_END_ANGLE) * i / SSHAPE_NDIV; // reverse
    cp[*node].x = center.x + SSHAPE_SCALE * ((1.0 + SSHAPE_WIDTH) * SSHAPE_XRAD * cos(theta));
    cp[*node].y = center.y + SSHAPE_SCALE * ((1.0 + SSHAPE_WIDTH) * SSHAPE_YRAD * sin(theta) + SSHAPE_LATTICE / 4.0);
    (*node) ++;
  }
  
  for (int i = 0; i < p; i ++) { // wrap p points to the first p ones
    double theta = SSHAPE_TOP_STA_ANGLE + (SSHAPE_TOP_END_ANGLE - SSHAPE_TOP_STA_ANGLE) * i / SSHAPE_NDIV;
    cp[*node].x = center.x + SSHAPE_SCALE * ((1.0 - SSHAPE_WIDTH) * SSHAPE_XRAD * cos(theta));
    cp[*node].y = center.y + SSHAPE_SCALE * ((1.0 - SSHAPE_WIDTH) * SSHAPE_YRAD * sin(theta) + SSHAPE_LATTICE / 4.0);
    (*node) ++;
  }

}

void isolap2d_input_sshape(Isolap2d *p)
{
  const int degree = SSHAPE_SPLINE_DEGREE;
  const int num_per_sshape = 4 * SSHAPE_NDIV + degree; // number of CP per curve

  const int nx = SSHAPE_NUM_FOR_X;
  const int ny = SSHAPE_NUM_FOR_Y;

  p->ncc = nx * ny; // number of closed curves
  INFO("ncc=%d\n", p->ncc);

  p->cc = closedcurve_allocate(p->ncc);

  int icc = 0; // initialise the index for closed curve

  for (int j = 0; j < ny; j ++) {
    for (int i = 0; i < nx; i ++) {

      (p->cc)[icc].p = degree;
      (p->cc)[icc].n = num_per_sshape;

      (p->cc)[icc].cp = closedcurve_cp_alloc((p->cc)[icc].n); // control points

      double2 center;
      center.x = SSHAPE_LATTICE * i;
      center.y = SSHAPE_LATTICE * j;

      (p->cc)[icc].x = center; // store the center for the purpose of meshing etc

      int node = 0;

      sshape((p->cc)[icc].p, center, &node, (p->cc)[icc].cp);
 
      icc ++;
    }
  }
  
}

void isolap2d_input_sshape_info()
{
  INFO("SSHAPE_SPLINE_DEGREE=%d\n", SSHAPE_SPLINE_DEGREE);
  INFO("SSHAPE_LATTICE=%f\n", SSHAPE_LATTICE);
  INFO("SSHAPE_SCALE=%f\n", SSHAPE_SCALE);
  INFO("SSHAPE_WIDTH=%f\n", SSHAPE_WIDTH);
  INFO("SSHAPE_NDIV=%d\n", SSHAPE_NDIV);
  INFO("SSHAPE_XRAD=%f\n", SSHAPE_XRAD);
  INFO("SSHAPE_YRAD=%f\n", SSHAPE_YRAD);

  INFO("SSHAPE_NUM_FOR_X=%d\n", SSHAPE_NUM_FOR_X);
  INFO("SSHAPE_NUM_FOR_Y=%d\n", SSHAPE_NUM_FOR_Y);
}
