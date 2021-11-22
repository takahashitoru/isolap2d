#include "quadtree.h"

int imleaf(const Cell c)
{
  if (c.child[0] == NULL_CELL && c.child[1] == NULL_CELL && c.child[2] == NULL_CELL && c.child[3] == NULL_CELL) {
    return QUADTREE_TRUE;  // leaf
  } else {
    return QUADTREE_FALSE; // non-leaf
  }
}

/*
  ibrnch=(XY)_2
  XDIR(0)=-1 YDIR(0)=-1
  XDIR(1)=-1 YDIR(1)= 1
  XDIR(2)= 1 YDIR(2)=-1
  XDIR(3)= 1 YDIR(3)= 1
*/

#define XDIR(i) ( 2 * (((i) >> 1) & 1) - 1 )
#define YDIR(i) ( 2 * (((i)     ) & 1) - 1 )

static void chidir(const int ibrnch, double2 *dir)
{
  //////////////////////////////////
  ASSERT(0 <= ibrnch && ibrnch < 4);
  //////////////////////////////////
  dir->x = XDIR(ibrnch);
  dir->y = YDIR(ibrnch);
}


static int ibrnch(const double2 c, const double2 x)
{
  int s = 0;

  if (x.y <= c.y) {
    s += 0;
  } else { // c.y < x.y
    s += 1;
  }

  if (x.x <= c.x) {
    s += 0;
  } else { // c.x < x.x
    s += 2;
  }

  return s;
}


static int this_is_a_leaf(const int minlev, const int maxepc, const int level, const Cell c)
{
  if (c.btail - c.bhead + 1 <= maxepc
      && c.ctail - c.chead + 1 <= maxepc
      && level >= minlev) {
    return QUADTREE_TRUE; // this is a leaf
  } else {
    return QUADTREE_FALSE;  // this is not a leaf
  }
}


#if(1) // sort without sort

static void sort_list(const double2 center, const int head, const int tail, const double2 *x, int *array, int *counts, int *heads)
{
  const int num = tail - head + 1;

  int *lists = (int *)malloc(4 * num * sizeof(int)); // lists[4][num]

  for (int s = 0; s < 4; s ++) {
    counts[s] = 0; // initialise
  } 

  for (int i = head; i <= tail; i ++) {
    const int s = ibrnch(center, x[array[i]]); // obtain the sibling index; 0<=s<4
    lists[s * num + counts[s]] = array[i];
    counts[s] ++;
  }

  int p = head;
  for (int s = 0; s < 4; s ++) {
    heads[s] = p; // heads[0]=head, heads[1]=head+counts[0], ...
    for (int n = 0; n < counts[s]; n ++) {
      array[p] = lists[s * num + n];
      p ++;
    }
  }

  free(lists);

}

#else

typedef struct sort{
  int e;
  double2 c;
  double2 x;
} Sort;


static int comp(const void *pa, const void *pb)
{
  Sort a = *(Sort *)pa;
  Sort b = *(Sort *)pb;

  /* b.c=a.c */

  int ia = ibrnch(a.c, a.x);
  int ib = ibrnch(b.c, b.x);

  if (ia < ib) {
    return - 1;
  } else if (ia > ib) {
    return 1;
  } else {
    return 0;
  }
}


static void sort_list(const double2 center, const int head, const int tail, const double2 *x, int *array)
{
#ifdef DISABLE_QSORT
  for (int i = head; i <= tail; i ++) {
    for (int j = i + 1; j <= tail; j ++) {
      if (ibrnch(center, x[array[i]]) > ibrnch(center, x[array[j]])) {
	int itmp = array[i];
	array[i] = array[j];
	array[j] = itmp;
      }
    }
  }
#else
  int n = tail - head + 1; // number of points in this cell (n>=0)
  Sort *iwork = (Sort *)malloc(n * sizeof(Sort));

  n = 0; // reset counter
  for (int i = head; i <= tail; i ++) {
    iwork[n].e = array[i];
    iwork[n].c = center;
    iwork[n].x = x[array[i]];
    n ++;
  }
  qsort(iwork, n, sizeof(Sort), &comp);
  n = 0; // reset counter
  for (int i = head; i <= tail; i ++) {
    array[i] = iwork[n].e;
    n ++;
  }

  free(iwork);
#endif
}

#endif


void quadtree(const int maxdep, const int maxcel, const int maxepc, Cell *c,
	      int *maxlev, const int minlev, int *ncell, int *levsta, int *levend, 
	      double *celeng, int *barray, int *carray,
	      const double2 rotcnt, const double rotlen,
	      //120121	      const int p, const int n, const double2 *xb, const double2 *xc)
	      const int nbasis, const int ncol, const double2 *xb, const double2 *xc)
{
  /* Initialise the array for basis functions and the array for
     collocation points */

  //120121  for (int i = 0; i < n; i ++) {
  for (int i = 0; i < nbasis; i ++) {
    barray[i] = i;
  }

  //120121  for (int i = 0; i < n - p; i ++) {
  for (int i = 0; i < ncol; i ++) {
    carray[i] = i;
  }

  /* Set up the root cell */

  *ncell = ROOT_CELL;
  c[*ncell].bhead = 0;
  //120121  c[*ncell].btail = n - 1;
  c[*ncell].btail = nbasis - 1;
  c[*ncell].chead = 0;
  //120121  c[*ncell].ctail = n - p - 1;
  c[*ncell].ctail = ncol - 1;
  c[*ncell].parent = NULL_CELL;
  c[*ncell].center.x = rotcnt.x;
  c[*ncell].center.y = rotcnt.y;
  (*ncell) ++; // update

  /* Set up level of 0 */

  celeng[0] = rotlen;
  levsta[0] = 0;
  levend[0] = 0;
     
  /* Create cells */

  for (int level = 0; level < maxdep; level ++) { // up to level maxdep-1

    /* Initialise the flag to determine to go down to the next level */

    int iflag = QUADTREE_FALSE;

    /* Loop over cells in this level */

    for (int icell = levsta[level]; icell <= levend[level]; icell ++) {
      
      /* Check if this cell is a leaf or not */
      
      if (this_is_a_leaf(minlev, maxepc, level, c[icell]) == QUADTREE_FALSE) {
	
	/* If there is at least one non-leaf in this level, the flag
	   is turned on */

	iflag = QUADTREE_TRUE;

#if(1) // sort without sort
	/* Sort the bases and collocation lists for this cell by
	   sibling index */

	int icountb[4], koheadb[4];
	int icountc[4], koheadc[4];

	sort_list(c[icell].center, c[icell].bhead, c[icell].btail, xb, barray, icountb, koheadb);
	sort_list(c[icell].center, c[icell].chead, c[icell].ctail, xc, carray, icountc, koheadc);

#else

	/* Sort the bases and collocation lists for this cell by
	   sibling index */

	sort_list(c[icell].center, c[icell].bhead, c[icell].btail, xb, barray);
	sort_list(c[icell].center, c[icell].chead, c[icell].ctail, xc, carray);

	/* Intialise the numbers of bases and collocation points in
	   the i-th sibling */

	int icountb[4], koheadb[4];
	int icountc[4], koheadc[4];

	for (int i = 0; i < 4; i ++) {
	  icountb[i] = 0;
	  icountc[i] = 0;
	  koheadb[i] = c[icell].bhead;
	  koheadc[i] = c[icell].chead;
	}
	
	for (int i = c[icell].btail; i >= c[icell].bhead; i --) {
	  const int itmp = ibrnch(c[icell].center, xb[barray[i]]);
	  icountb[itmp] ++;
	  koheadb[itmp] = i;
	}

	for (int i = c[icell].ctail; i >= c[icell].chead; i --) {
	  const int itmp = ibrnch(c[icell].center, xc[carray[i]]);
	  icountc[itmp] ++;
	  koheadc[itmp] = i;
	}
#endif

	/* Create children of this cell */

	for (int i = 0; i < 4; i ++) {

	  if (icountb[i] > 0 || icountc[i] > 0) {
	    
	    c[*ncell].parent = icell;
	    c[icell].child[i] = *ncell;

	    c[*ncell].bhead = koheadb[i];
	    c[*ncell].btail = koheadb[i] + icountb[i] - 1;

	    c[*ncell].chead = koheadc[i];
	    c[*ncell].ctail = koheadc[i] + icountc[i] - 1;

	    double2 dir;
	    chidir(i, &dir);
	    c[*ncell].center.x = c[icell].center.x + dir.x * celeng[level] / 4;
	    c[*ncell].center.y = c[icell].center.y + dir.y * celeng[level] / 4;

	    (*ncell) ++;

	    if (*ncell > maxcel) {
	      MESG("maxcel is too small. Exit.\n");
	      exit(EXIT_FAILURE);
	    }

	  } else { // icountb[i] == 0 && icountc[i] == 0

	    c[icell].child[i] = NULL_CELL;

	  }
	}

      } else { // icell is a leaf

	for (int i = 0; i < 4; i ++) {
	  c[icell].child[i] = NULL_CELL;
	}

      }
    } // icell
    

    if (iflag == QUADTREE_FALSE) {

      /* Since all the cells in this level are leaves, we can regard
	 this level as the maximum level and can exit */

      *maxlev = level;

      return;

    } else {

      /* Set up for the next level */

      celeng[level + 1] = celeng[level] / 2;
      levsta[level + 1] = levend[level] + 1;
      levend[level + 1] = *ncell - 1; // the last cell index is not *ncell

    }

  } // level

  MESG("maxdep is too small. Exit.\n");
  exit(EXIT_FAILURE);

}
