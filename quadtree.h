#ifndef QUADTREE_H
#define QUADTREE_H

#include <stdio.h>
#include <stdlib.h>

#include "debugmacros.h"
#include "double2.h"
#include "cell.h"

/* Constants */

#define QUADTREE_TRUE  (1)
#define QUADTREE_FALSE (0)

#define NULL_CELL (- 1)
#define ROOT_CELL (0)

#define NULL_LEVEL (- 1)

/* Structure */

typedef struct {
  int maxlev;
  int minlev;
  int ncell;
  int *barray; // list of basis functions
  int *carray; // list of collocation points;
  int *levsta;
  int *levend;
  double *celeng;
  Cell *cell;
} Tree;

/* Function prototypes */

#ifdef __cplusplus
extern "C" {
#endif
  int imleaf(const Cell c);
  void quadtree(const int maxdep, const int maxcel, const int maxepc, Cell *c,
		int *maxlev, const int minlev, int *ncell, int *levsta, int *levend, 
		double *celeng, int *barray, int *carray,
		const double2 rotcnt, const double rotlen,
		//120121		const int p, const int n, const double2 *xb, const double2 *xc);
		const int nbase, const int ncol, const double2 *xb, const double2 *xc);
#ifdef __cplusplus
}
#endif

#endif /* QUADTREE_H */
