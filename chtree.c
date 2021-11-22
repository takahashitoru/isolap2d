#include "isolap2d.h"

#define CHTREE_TMPFILE "CHTREE_cell.txt"
#define CHTREE_TMPFILE2 "CHTREE_cell.label"
#define CHTREE_TMPFILE3 "CHTREE_leaf.txt"
#define CHTREE_TMPFILE4 "CHTREE_leaf.label"

static void draw_a_cell_by_gnuplot(FILE *fp, int level, Cell c, double *celeng)
{
  double tmp = celeng[level] / 2;
  fprintf(fp, "%15.7e %15.7e\n", c.center.x + tmp, c.center.y + tmp);
  fprintf(fp, "%15.7e %15.7e\n", c.center.x - tmp, c.center.y + tmp);
  fprintf(fp, "%15.7e %15.7e\n", c.center.x - tmp, c.center.y - tmp);
  fprintf(fp, "%15.7e %15.7e\n", c.center.x + tmp, c.center.y - tmp);
  fprintf(fp, "%15.7e %15.7e\n", c.center.x + tmp, c.center.y + tmp);
  fprintf(fp, "\n"); /* brank line */
}

static void print_a_cell_index_by_gnuplot(FILE *fp, Cell *c, int icell)
{
  fprintf(fp, "set label \"%d\" at %e, %e center\n", icell, c[icell].center.x, c[icell].center.y);
}

static void chtree(int *levsta, int *levend, int maxlev, Cell *c, double *celeng)
{
  FILE *fp, *fp2, *fp3, *fp4;

  if ((fp = fopen(CHTREE_TMPFILE, "w")) == NULL) {
    INFO("Fail to open. Exit.\n");
    exit(EXIT_FAILURE);
  }

  if ((fp2 = fopen(CHTREE_TMPFILE2, "w")) == NULL) {
    INFO("Fail to open. Exit.\n");
    exit(EXIT_FAILURE);
  }

  if ((fp3 = fopen(CHTREE_TMPFILE3, "w")) == NULL) {
    INFO("Fail to open. Exit.\n");
    exit(EXIT_FAILURE);
  }

  if ((fp4 = fopen(CHTREE_TMPFILE4, "w")) == NULL) {
    INFO("Fail to open. Exit.\n");
    exit(EXIT_FAILURE);
  }

  for (int level = 0; level <= maxlev; level++) {
    for (int icell = levsta[level]; icell <= levend[level]; icell++) {

      draw_a_cell_by_gnuplot(fp, level, c[icell], celeng);
      print_a_cell_index_by_gnuplot(fp2, c, icell);

      if (imleaf(c[icell]) == QUADTREE_TRUE) {
	draw_a_cell_by_gnuplot(fp3, level, c[icell], celeng);
	print_a_cell_index_by_gnuplot(fp4, c, icell);
      }

    }
  }     

  fclose(fp);
  fclose(fp2);
  fclose(fp3);
  fclose(fp4);
}

void isolap2d_chtree(Isolap2d *p)
{
  Tree *tree = p->tree;
  chtree(tree->levsta, tree->levend, tree->maxlev, tree->cell, tree->celeng);
}
