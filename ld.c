#include "isolap2d.h"

#if !defined(LOAD_PARAMS_FROM_FILE)

#ifndef NGAUSS
#define NGAUSS 4
#endif
#ifndef NINC
#define NINC 300
#endif
#ifndef EPS
#define EPS 1.0e-12 //1.0e-12
#endif
#define MAXL 1000
#define JPRE 0
#ifndef ITMAX
#define ITMAX 10
#endif
#ifndef TOL
#define TOL 1.0e-5
#endif
#ifndef INITX
#define INITX 0
#endif
#if defined(FMM)
#define RATIO 1.01
#define MAXDEP 20
#define MINDEP 2
#define MAXCEL 1000000
#ifndef MAXEPC
#define MAXEPC 30
#endif
#ifndef NBETA
#define NBETA 1
#endif
#ifndef NTERM
#define NTERM 10
#endif
#endif

void isolap2d_ld(Isolap2d *p)
{
  p->params = (Params *)malloc(sizeof(Params));
  Params *pa = p->params; // alias

  pa->ngauss = NGAUSS; // # of Gaussian quadrature points per element [Unused]
  pa->ninc = NINC;     // ID of the external field (ignore -> 0)
  pa->eps = EPS;       // small number
  pa->maxl = MAXL;     // Dimension of Krylov subspace
  pa->jpre = JPRE;     // ID to choose preconditioning (no prec. -> 0)
  pa->itmax = ITMAX;   // Maximum # of iterations
  pa->tol = TOL;       // Torlerance
  pa->initx = INITX;   // Set 0 if you give the null initial guess.
#if defined(FMM)
  pa->ratio = RATIO;   // FMM
  pa->maxdep = MAXDEP; // FMM
  pa->mindep = MINDEP; // FMM
  pa->maxcel = MAXCEL; // FMM
  pa->maxepc = MAXEPC; // FMM; maximum number of collocation points or bases per cell
  pa->nbeta = NBETA;   // FMM; order of neighbors
  pa->nterm = NTERM;   // FMM; number of multipoles
#endif

}

#else

void isolap2d_ld(char *filename, Isolap2d *p)
{
  FILE *fp = fopen(filename, "r");
  
  if (fp == NULL) {
    INFO("Fail to open %s; exit.\n", filename);
    exit(1);
  }

  p->params = (Params *)malloc(sizeof(Params));
  Params *pa = p->params; // alias

  fscanf(fp, "%d", &(pa->ngauss)); // # of Gaussian quadrature points per element [Unused]
  fscanf(fp, "%d", &(pa->ninc));   // ID of the external field (ignore -> 0)
  fscanf(fp, "%lf", &(pa->eps));   // small number
  fscanf(fp, "%d", &(pa->maxl));   // Dimension of Krylov subspace
  fscanf(fp, "%d", &(pa->jpre));   // ID to choose preconditioning (no prec. -> 0)
  fscanf(fp, "%d", &(pa->itmax));  // Maximum # of iterations
  fscanf(fp, "%lf", &(pa->tol));   // Torlerance
  fscanf(fp, "%d", &(pa->initx));  // Set 0 if you give the null initial guess.
#if defined(FMM)
  fscanf(fp, "%lf", &(pa->ratio));
  fscanf(fp, "%d", &(pa->maxdep));
  fscanf(fp, "%d", &(pa->mindep));
  fscanf(fp, "%d", &(pa->maxcel));
  fscanf(fp, "%d", &(pa->maxepc));
  fscanf(fp, "%d", &(pa->nbeta));
  fscanf(fp, "%d", &(pa->nterm));
#endif

  fclose(fp);
}

#endif

void isolap2d_ld_check(const Isolap2d *p)
{
  const Params *pa = p->params;

  INFO("ngauss=%d\n", pa->ngauss);
  INFO("ninc=%d\n", pa->ninc);
  INFO("eps=%e\n", pa->eps);
  INFO("maxl=%d\n", pa->maxl);
  INFO("jpre=%d\n", pa->jpre);
  INFO("itmax=%d\n", pa->itmax);
  INFO("tol=%e\n", pa->tol);
  INFO("initx=%d\n", pa->initx);
#if defined(FMM)
  INFO("ratio=%e\n", pa->ratio);
  INFO("maxdep=%d\n", pa->maxdep);
  INFO("mindep=%d\n", pa->mindep);
  INFO("maxcel=%d\n", pa->maxcel);
  INFO("maxepc=%d\n", pa->maxepc);
  INFO("nbeta=%d\n", pa->nbeta);
  INFO("nterm=%d\n", pa->nterm);
#endif
}

void isolap2d_ld_free(Isolap2d *p)
{
  Params *pa = p->params;
  free(pa);
}
