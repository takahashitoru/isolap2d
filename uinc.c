#include "isolap2d.h"

static double uinc(const int ninc, const double2 x)
{
  double u;

  switch (ninc) {
  case 100:
    u = x.x;
    break;

  case 200:
    u = x.y;
    break;

  case 300:
    u = x.y;
    break;

    /* Regist incidental wave field */

    //  case 300:
    //    u = ...;
    //    break;
    
  default:
    MESG("Unregistered ninc was used. Exit.\n");
    exit(EXIT_FAILURE);
  }
  
  return u;
}


void mkuinc(const int ncol, const int ninc, const double2 *xc, double *bvec)
{
  if (ninc != NULL_UINC) {
    for (int ix = 0; ix < ncol; ix ++) {
      const double u = uinc(ninc, xc[ix]);
      bvec[ix] += u;
    }
  }
}

void isolap2d_mkuinc(const Isolap2d *p, const int ncol, double *bvec)
{
  const Params *pa = p->params;
  mkuinc(ncol, pa->ninc, p->xc, bvec);
}
