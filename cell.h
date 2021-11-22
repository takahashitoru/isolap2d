#ifndef CELL_H
#define CELL_H

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include "double2.h"
#include "dcomplex.h"

typedef struct {
  //  int head, head2;
  //  int tail, tail2;
  int bhead, btail; // for basis functions
  int chead, ctail; // for collocation points
  int parent;
  int child[4];
  double2 center;
  int nneigh;
  int *neighbors;
  int ninter;
  int *interactions;
  //  dcomplex *mcoe[2];
  //  dcomplex *lcoe[2];
  //  dcomplex *fcoe[2];
  //  dcomplex *hcoe[2];
  dcomplex *mcoe;
  dcomplex *lcoe;
} Cell;

///* Macro to count the number of boundary elements in a cell */
//#define GET_NUMBER_OF_ELEMENTS(icell) (c[icell].tail - c[icell].head + 1)
//
///* Macro to judge the cell has no boundary elements */
//#define IS_EMPTY(icell) (GET_NUMBER_OF_ELEMENTS(icell) <= 0)
//#define IS_NOT_EMPTY(icell) (!IS_EMPTY(icell))
//
//#define IS_LEAF(level, icell) (imleaf(c[icell]))
//#define IS_APPROPRIATE(level, icell) (TRUE) // always true

#endif /* CELL_H */

