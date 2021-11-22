#ifndef ISOLAP2D_DEFS_H
#define ISOLAP2D_DEFS_H

#include "double2.h"  // defines double2
#include "gaussone.h" // defines GaussOne
#include "quadtree.h" // defines Tree
#include "closedcurve.h" // defines Closedcurve


typedef struct {
  int ngauss;
  int ninc;
  double eps;
  int maxl;
  int jpre;
  int itmax;
  double tol;
  int initx;

  double ratio; // FMM
  int maxdep; // FMM
  int mindep; // FMM
  int maxcel; // FMM
  int maxepc; // FMM
  int nbeta; // FMM
  int nterm; // FMM

} Params;

typedef struct {
  //120121  int n;            // number of control points
  //120121  int p;            // order of base functions; two or more.
  //120121  double *t;        // knots
  //120121  double2 *cp;      // control points
  int ncc;          // number of closed curves
  Closedcurve *cc; // closed curves
  double *tc;       // knots corresponding to collocation points
  double2 *xc;      // collocation points
  int *ccc;         // index of closed curve of collocation points
  int *xoff;        // offset to xvec
  double2 *xb;      // centre of basis functions
  int *bcc;         // index of closed curve of basis
  int *bj;          // index of local basis
  double *freeterm; // free term
  GaussOne *gauss;  // one-dimensional Gausssian integral formula
  Params *params;
  Tree *tree;
} Isolap2d;

#endif /* ISOLAP2D_DEFS_H */
