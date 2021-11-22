#include "isolap2d.h"

#define BIG_NUMBER (1.0e+200)

void isolap2d_setup_root(const int nbasis, const int ncol, const double ratio, const double2 *xb, const double2 *xc, double2 *rotcnt, double *rotlen)
{

  double2 dmax, dmin;
  dmax.x = dmax.y = - BIG_NUMBER;
  dmin.x = dmin.y = BIG_NUMBER;

  for (int i = 0; i < nbasis; i ++) { // base functions
    dmax.x = fmax(dmax.x, xb[i].x);
    dmax.y = fmax(dmax.y, xb[i].y);
    dmin.x = fmin(dmin.x, xb[i].x);
    dmin.y = fmin(dmin.y, xb[i].y);
  }

  for (int i = 0; i < ncol; i ++) { // collocation points
    dmax.x = fmax(dmax.x, xc[i].x);
    dmax.y = fmax(dmax.y, xc[i].y);
    dmin.x = fmin(dmin.x, xc[i].x);
    dmin.y = fmin(dmin.y, xc[i].y);
  }
      
  (*rotcnt).x = (dmax.x + dmin.x) / 2.0;
  (*rotcnt).y = (dmax.y + dmin.y) / 2.0;

  double2 tmp;
  tmp.x = dmax.x - dmin.x;
  tmp.y = dmax.y - dmin.y;
  
  *rotlen = fmax(tmp.x, tmp.y) * ratio;

  INFO("rotcnt=%e %e\n", (*rotcnt).x, (*rotcnt).y);
  INFO("rotlen=%e\n", *rotlen);
}


#ifdef __ICC
#pragma intel optimization_level 0
#endif
//120124void isolap2d_make_tree(Isolap2d *p, const double2 rotcnt, const double rotlen, const int nbasis, const int ncol)
void isolap2d_make_tree(Isolap2d *p, const double2 rotcnt, const double rotlen,
			const int nbasis, const int ncol, const double2 *xb, const double2 *xc)
{
  Params *pa = p->params;

  p->tree = (Tree *)malloc(sizeof(Tree)); // allocate

  Tree *tree = p->tree;

  tree->barray = (int *)malloc(nbasis * sizeof(int));                   // barray[0:nbasis-1]
  tree->carray = (int *)malloc(ncol * sizeof(int));                     // carray[0:ncol-1]
  tree->levsta = (int *)malloc(((pa->maxdep) + 1) * sizeof(int));       // levsta[0:maxdep]
  tree->levend = (int *)malloc(((pa->maxdep) + 1) * sizeof(int));       // levend[0:maxdep]
  tree->celeng = (double *)malloc(((pa->maxdep) + 1) * sizeof(double)); // celeng[0:maxdep]
  tree->cell = (Cell *)malloc((pa->maxcel) * sizeof(Cell));             // cell[0:maxcel-1]

  if (pa->nbeta == 1 || pa->nbeta == 2) {
    tree->minlev = 2;
  } else {
    INFO("Unsupposed nbeta=%d\n. Exit.\n", pa->nbeta);
    exit(EXIT_FAILURE);
  }

  INFO("minlev=%d\n", tree->minlev);

  quadtree(pa->maxdep, pa->maxcel, pa->maxepc, tree->cell,
	   &(tree->maxlev), tree->minlev, &(tree->ncell), tree->levsta, tree->levend,
	   tree->celeng, tree->barray, tree->carray,
	   rotcnt, rotlen,
	   //120124	   nbasis, ncol, p->xb, p->xc);
	   nbasis, ncol, xb, xc);

  INFO("maxlev=%d\n", tree->maxlev);
  INFO("ncell=%d\n", tree->ncell);

}

void isolap2d_free_tree(Isolap2d *p)
{
  Tree *tree = p->tree;
  free(tree->barray);
  free(tree->carray);
  free(tree->levsta);
  free(tree->levend);
  free(tree->celeng);
  free(tree->cell);
  free(tree);
}
