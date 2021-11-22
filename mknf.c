#include "isolap2d.h"

static void mknf(const int nbeta, const int maxlev, const int *levsta, const int *levend, const double *celeng, Cell *c)
{
  const int root = ROOT_CELL;
  c[root].nneigh = 1;
  c[root].neighbors[0] = root;
  c[root].ninter = 0;

  for (int level = 1; level <= maxlev; level ++) {

    double boui = (nbeta + 0.5) * celeng[level];
    boui *= boui;

    //////////////////////////////////////////////////////////////////////////
    //    DBG("level=%d # of cells=%d\n", level, levend[level] - levsta[level] + 1);
    //////////////////////////////////////////////////////////////////////////

#if defined(_OPENMP)
#pragma omp parallel for // OpenMP DEFINED LOOP WAS PARALLELIZED.
#endif
    for (int icell = levsta[level]; icell <= levend[level]; icell ++) {
      c[icell].nneigh = 0; // initialise the number of near neighbour cells
      c[icell].ninter = 0; // initialise the number of interacting cells

      const int iparent = c[icell].parent;
	
      for (int j = 0; j < c[iparent].nneigh; j ++) {
	const int iuncle = c[iparent].neighbors[j];

	for (int i = 0; i < 4; i ++) {
	  const int icousn = c[iuncle].child[i];

	  if (icousn != NULL_CELL) { // only if i-th child of iuncle exists

	    double2 tmp;
	    tmp = sub2(c[icell].center, c[icousn].center);
	    tmp.x *= tmp.x;
	    tmp.y *= tmp.y;
		
	    if ((tmp.x < boui) && (tmp.y < boui)) {
	      c[icell].neighbors[c[icell].nneigh] = icousn;
	      (c[icell].nneigh) ++;
	    } else {
	      c[icell].interactions[c[icell].ninter] = icousn;
	      (c[icell].ninter) ++;
	    }

	  }
	} // i
      } // j
    } // icell
  } // level
}


static void check_mknf(const int maxlev, const int *levsta, const int *levend, const Cell *c)
{
  for (int level = 0; level <= maxlev; level ++) {
    for (int icell = levsta[level]; icell <= levend[level]; icell ++) {
      DBG("level=%d icell=%d leaf=%d nneigh=%d ninter=%d\n",
      	  level, icell, imleaf(c[icell]), c[icell].nneigh, c[icell].ninter);
#if(1)
      MSG("interactions:");
      for (int iinter = 0; iinter < c[icell].ninter; iinter ++) {
	fprintf(stderr, " %d", c[icell].interactions[iinter]);	
      }
      fprintf(stderr, "\n");
      MSG("neighbors   :");
      for (int ineigh = 0; ineigh < c[icell].nneigh; ineigh ++) {
	fprintf(stderr, " %d", c[icell].neighbors[ineigh]);	
      }
      fprintf(stderr, "\n");
#endif
    }
  }
}


#if defined(__ICC)
#pragma intel optimization_level 0
#endif
void isolap2d_mknf(Isolap2d *p)
{
  Params *pa = p->params;
  Tree *tree = p->tree;

  /* Allocate neighbours and interactions of every cells */

  const int max_neighbors = SQUARE(2 * pa->nbeta + 1);
  const int max_interactions = SQUARE(4 * pa->nbeta + 2) - SQUARE(2 * pa->nbeta + 1);

  INFO("max_neighbors=%d\n", max_neighbors);
  INFO("max_interactions=%d\n", max_interactions);

  for (int level = 0; level <= tree->maxlev; level ++) {
    for (int icell = (tree->levsta)[level]; icell <= (tree->levend)[level]; icell ++) {
      (tree->cell)[icell].neighbors = (int *)malloc(max_neighbors * sizeof(int));
      (tree->cell)[icell].interactions = (int *)malloc(max_interactions * sizeof(int));
    }
  }

  mknf(pa->nbeta, tree->maxlev, tree->levsta, tree->levend, tree->celeng, tree->cell);
}


void isolap2d_mknf_check(const Isolap2d *p)
{
  Tree *tree = p->tree;
  check_mknf(tree->maxlev, tree->levsta, tree->levend, tree->cell);
}


void isolap2d_mknf_free(Isolap2d *p)
{
  Tree *tree = p->tree;

  for (int level = 0; level <= tree->maxlev; level ++) {
    for (int icell = (tree->levsta)[level]; icell <= (tree->levend)[level]; icell ++) {
      free((tree->cell)[icell].neighbors);
      free((tree->cell)[icell].interactions);
    }
  }
}
