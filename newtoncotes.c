#include "newtoncotes.h"

static void newtoncotes(const int n, double *w)
{
  /* See 25.4.13--25.4.20 in A&S */

  if (n == 4) {
    const double c = 3.0 / 8.0;
    w[0] = c * 1;
    w[1] = c * 3;
    w[2] = c * 3;
    w[3] = c * 1;

  } else if (n == 5) {
    const double c = 2.0 / 45.0;
    w[0] = c *  7;
    w[1] = c * 32;
    w[2] = c * 12;
    w[3] = c * 32;
    w[4] = c *  7;

  } else if (n == 6) {
    const double c = 5.0 / 288.0;
    w[0] = c * 19;
    w[1] = c * 75;
    w[2] = c * 50;
    w[3] = c * 50;
    w[4] = c * 75;
    w[5] = c * 19;

  } else if (n == 7) {
    const double c = 1.0 / 140.0;
    w[0] = c *  41;
    w[1] = c * 216;
    w[2] = c *  27;
    w[3] = c * 272;
    w[4] = c *  27;
    w[5] = c * 216;
    w[6] = c *  41;

  } else if (n == 8) {
    const double c = 7.0 / 17280.0;
    w[0] = c *  751;
    w[1] = c * 3577;
    w[2] = c * 1323;
    w[3] = c * 2989;
    w[4] = c * 2989;
    w[5] = c * 1323;
    w[6] = c * 3577;
    w[7] = c *  751;

  } else if (n == 9) {
    const double c = 4.0 / 14175.0;
    w[0] = c * 989;
    w[1] = c * 5888;
    w[2] = c * (- 928);
    w[3] = c * 10496;
    w[4] = c * (- 4540);
    w[5] = c * 10496;
    w[6] = c * (- 928);
    w[7] = c * 5888;
    w[8] = c * 989;

  } else if (n == 10) {
    const double c = 9.0 / 89600.0;
    w[0] = c *  2857;
    w[1] = c * 15741;
    w[2] = c *  1080;
    w[3] = c * 19344;
    w[4] = c *  5778;
    w[5] = c *  5778;
    w[6] = c * 19344;
    w[7] = c *  1080;
    w[8] = c * 15741;
    w[9] = c *  2857;

  } else if (n == 11) {
    const double c = 5.0 / 299376.0;
    w[ 0] = c *  16067;
    w[ 1] = c * 106300;
    w[ 2] = c * (- 48525);
    w[ 3] = c * 272400;
    w[ 4] = c * (- 260550);
    w[ 5] = c * 427368;
    w[ 6] = c * (- 260550);
    w[ 7] = c * 272400;
    w[ 8] = c * (- 48525);
    w[ 9] = c * 106300;
    w[10] = c *  16067;

  } else {
    
    INFO("Invalid parameter: n=%d\n. Exit.", n);
    exit(EXIT_FAILURE);

  }

}


void newtoncotes_load(const int n, Newtoncotes **p)
{
  *p = (Newtoncotes *)malloc(sizeof(Newtoncotes));
  (*p)->n = n;
  (*p)->h = 1.0 / (n - 1);
  (*p)->w = (double *)malloc(n * sizeof(double));

  newtoncotes(n, (*p)->w);
}

void newtoncotes_check(const Newtoncotes *p)
{
  INFO("n=%d h=%le\n", p->n, p->h);
  for (int i = 0; i < p->n; i ++) {
    INFO("i=%d (x=%le) w=%le\n", i, p->h * i, (p->w)[i]);
  }
}

void newtoncotes_free(Newtoncotes *p)
{
  free(p->w);
}
