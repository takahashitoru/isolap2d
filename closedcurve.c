#include "closedcurve.h"

Closedcurve *closedcurve_allocate(const int ncc)
{
  return (Closedcurve *)malloc(ncc * sizeof(Closedcurve));
}


void closedcurve_free(Closedcurve *cc)
{
  free(cc);
}


double2 *closedcurve_cp_alloc(const int n)
{
  /* Allocate n control points */
  
  if (n <= 0) {
    fprintf(stderr, "Invaild n=%d. Exit.\n", n);
    exit(EXIT_FAILURE);
  }

  return (double2 *)malloc(n * sizeof(double2));
}


void closedcurve_cp_free(double2 *cp)
{
  free(cp);
}


double *closedcurve_knot_alloc(const int p, const int n)
{
  /* Allocate n+p+1 knots for a closed curve; t[0:n+p] */
  
  if (p < 2) {
    fprintf(stderr, "Invaild p=%d. Exit.\n", p);
    exit(EXIT_FAILURE);
  }
  
  if (n <= 0) {
    fprintf(stderr, "Invaild n=%d. Exit.\n", n);
    exit(EXIT_FAILURE);
  }

  return (double *)malloc((n + p + 1) * sizeof(double));
}


void closedcurve_knot_free(double *t)
{
  free(t);
}


void closedcurve_knot_comp(const int p, const int n, double *t)
{
  /* Compute uniform knot values; t[0:n+p] */
  
  const double dt = 1.0 / (n + p); // uniform

  for (int i = 0; i < n + p + 1; i ++) {
    t[i] = i * dt;
  }

}


void closedcurve_knot_check(const int ncc, const Closedcurve *cc)
{
  FILE *fp = fopen("knot.txt", "w");
  FILE *fp2 = fopen("knot.label", "w");
  for (int icc = 0; icc < ncc; icc ++) {
    const int p = cc[icc].p;
    const int n = cc[icc].n;
    const double *t = cc[icc].t;
    const double2 *cp = cc[icc].cp;
    for (int i = p; i < n; i ++) { // only the boundary
      double2 c = closedcurve_position(p, n, t, cp, t[i]);
      fprintf(fp, "%f %f 0.0\n", c.x, c.y);
      fprintf(fp2, "set label 'T%d' at %f, %f, %f \n", i, c.x, c.y, 0.1 * i);
    }
    fprintf(fp, "\n\n");
  }
  fclose(fp);
  fclose(fp2);
  MSG("Created knot.txt and knot.label\n");
}


void closedcurve_draw_boundary(const int ncc, const Closedcurve *cc)
{
  FILE *fp = fopen("boundary.txt", "w");

  for (int icc = 0; icc < ncc; icc ++) {
    const int p = cc[icc].p;
    const int n = cc[icc].n;
    const double *t = cc[icc].t;
    const double2 *cp = cc[icc].cp;
    DBG("Domain of the curve [t_%d,t_%d)=[%f,%f)\n", p, n, t[p], t[n]);

    for (double x = t[p]; x < t[n]; x += (t[1] - t[0]) / 50) {
      const double2 c = closedcurve_position(p, n, t, cp, x);
      fprintf(fp, "%f %f 0.0\n", c.x, c.y);
    }  
    fprintf(fp, "\n\n");
  }
  fclose(fp);

  MSG("Created boundary.txt\n");
}

/* B-spline for p=2 */
#define N20(x) ( (x) * (x) / 2.0 )
#define N21(x) ( (- 2.0 * (x) * (x) + 6.0 * (x) - 3.0) / 2.0 )
#define N22(x) ( (3.0 - (x)) * (3.0 - (x)) / 2.0 )

#define DN20(x) ( (x) )
#define DN21(x) ( - 2.0 * (x) + 3.0 )
#define DN22(x) ( - 3.0 + (x) )

#define DDN20(x) ( 1.0 )
#define DDN21(x) ( - 2.0 )
#define DDN22(x) ( 1.0 )

#define DDDN20(x) ( 0.0 )
#define DDDN21(x) ( 0.0 )
#define DDDN22(x) ( 0.0 )

/* B-spline for p=3 */
#define N30(x) ( (x) * (x) * (x) / 6.0 )
#define N31(x) ( (- 3.0 * (x) * (x) * (x) + 12.0 * (x) * (x) - 12.0 * (x) + 4.0) / 6.0 )
#define N32(x) ( (3.0 * (x) * (x) * (x) - 24.0 * (x) * (x) + 60.0 * (x) - 44.0) / 6.0 )
#define N33(x) ( - ((x) - 4.0) * ((x) - 4.0) * ((x) - 4.0) / 6.0 )

#define DN30(x) ( (x) * (x) / 2.0 )
#define DN31(x) ( (- 3.0 * (x) * (x) + 8.0 * (x) - 4.0) / 2.0 )
#define DN32(x) ( (3.0 * (x) * (x) - 16.0 * (x) + 20.0) / 2.0 )
#define DN33(x) ( - ((x) - 4.0) * ((x) - 4.0) / 2.0 )

#define DDN30(x) ( (x) )
#define DDN31(x) ( - 3.0 * (x) + 4.0 )
#define DDN32(x) ( 3.0 * (x) - 8.0 )
#define DDN33(x) ( - (x) + 4.0 )

#define DDDN30(x) ( 1.0 )
#define DDDN31(x) ( - 3.0 )
#define DDDN32(x) ( 3.0 )
#define DDDN33(x) ( - 1.0 )

/* B-spline for p=4 (thanks to Y. Nakai) */
#define N40(t) ( (t) * (t) * (t) * (t) / 24.0 )
#define N41(t) ( (- 4.0 * (t) * (t) * (t) * (t) + 20.0 * (t) * (t) * (t) - 30.0 * (t) * (t) + 20.0 * (t) - 5.0) / 24.0 )
#define N42(t) ( (6.0 * (t) * (t) * (t) * (t) - 60.0 * (t) * (t) * (t) + 210.0 * (t) * (t) - 300.0 * (t) + 155.0) / 24.0 )
#define N43(t) ( (- 4.0 * (t) * (t) * (t) * (t) + 60.0 * (t) * (t) * (t) - 330.0 * (t) * (t) + 780.0 * (t) - 655.0) / 24.0 )
#define N44(t) ( ((t) - 5.0) * ((t) - 5.0) * ((t) - 5.0) * ((t) - 5.0) / 24.0 )

#define DN40(t) ( (t) * (t) * (t) / 6.0 )
#define DN41(t) ( (- 4.0 * (t) * (t) * (t) + 15.0 * (t) * (t) - 15.0 * (t) + 5.0) / 6.0 )
#define DN42(t) ( (6.0 * (t) * (t) * (t) - 45.0 * (t) * (t) + 105.0 * (t) - 75.0) / 6.0 )
#define DN43(t) ( (- 4.0 * (t) * (t) * (t) + 45.0 * (t) * (t) - 165.0 * (t) + 195.0) / 6.0 )
#define DN44(t) ( ((t) - 5.0) * ((t) - 5.0) * ((t) - 5.0) / 6.0 )

#define DDN40(t) ( (t) * (t) / 2.0 )
#define DDN41(t) ( (- 4.0 * (t) * (t) + 10.0 * (t) - 5.0) / 2.0 )
#define DDN42(t) ( (6.0 * (t) * (t) - 30.0 * (t) + 35.0) / 2.0 )
#define DDN43(t) ( (- 4.0 * (t) * (t) + 30.0 * (t) - 55.0) / 2.0 )
#define DDN44(t) ( ((t) - 5.0) * ((t) - 5.0) / 2.0 )

#define DDDN40(t) ( (t) )
#define DDDN41(t) ( - 4.0 * (t) + 5.0 )
#define DDDN42(t) ( 6.0 * (t) - 15.0 )
#define DDDN43(t) ( - 4.0 * (t) + 15.0 )
#define DDDN44(t) ( (t) - 5.0 )


double closedcurve_basis(const int p, const int n, const double *t, const int j, const double tnow)
{
  /* Compute b_{j,p}(t) for t=tnow */

  const double dt = 1.0 / (n + p); // interval for uniform knots
  const double d = (tnow - t[j]) / dt;

  double b;

  if (p == 2) {

    if (d < 0) {
      b = 0.0;
    } else if (d < 1) {
      b = N20(d);
    } else if (d < 2) {
      b = N21(d);
    } else if (d < 3) {
      b = N22(d);
    } else {
      b = 0.0;
    }

  } else if (p == 3) {

    if (d < 0) {
      b = 0.0;
    } else if (d < 1) {
      b = N30(d);
    } else if (d < 2) {
      b = N31(d);
    } else if (d < 3) {
      b = N32(d);
    } else if (d < 4) {
      b = N33(d);
    } else {
      b = 0.0;
    }

  } else if (p == 4) {

    if (d < 0) {
      b = 0.0;
    } else if (d < 1) {
      b = N40(d);
    } else if (d < 2) {
      b = N41(d);
    } else if (d < 3) {
      b = N42(d);
    } else if (d < 4) {
      b = N43(d);
    } else if (d < 5) {
      b = N44(d);
    } else {
      b = 0.0;
    }

  } else {
    INFO("p=%d is not implemented yet. Exit.\n", p);
    exit(EXIT_FAILURE);
  }

  return b;
}


double closedcurve_basis_derivative(const int p, const int n, const double *t, const int j, const double tnow)
{
  /* Compute the derivative of b_{j,p}(t) for t=tnow */

  const double dt = 1.0 / (n + p); // interval for uniform knots
  const double d = (tnow - t[j]) / dt;

  double db;

  if (p == 2) {

    if (d < 0) {
      db = 0.0;
    } else if (d < 1) {
      db = DN20(d);
    } else if (d < 2) {
      db = DN21(d);
    } else if (d < 3) {
      db = DN22(d);
    } else {
      db = 0.0;
    }

  } else if (p == 3) {

    if (d < 0) {
      db = 0.0;
    } else if (d < 1) {
      db = DN30(d);
    } else if (d < 2) {
      db = DN31(d);
    } else if (d < 3) {
      db = DN32(d);
    } else if (d < 4) {
      db = DN33(d);
    } else {
      db = 0.0;
    }

  } else if (p == 4) {

    if (d < 0) {
      db = 0.0;
    } else if (d < 1) {
      db = DN40(d);
    } else if (d < 2) {
      db = DN41(d);
    } else if (d < 3) {
      db = DN42(d);
    } else if (d < 4) {
      db = DN43(d);
    } else if (d < 5) {
      db = DN44(d);
    } else {
      db = 0.0;
    }

  } else {
    INFO("p=%d is not implemented yet. Exit.\n", p);
    exit(EXIT_FAILURE);
  }

  db /= dt; // due to the differentialtion of (tnow-t[j])/dt

  return db;
}


double closedcurve_potential(const int p, const int n, const double *t, const double *xvec, const double tnow)
{
  const double dt = 1.0 / (n + p); // interval for uniform knots
  int k = (int)(tnow / dt); // current point is in the interval [t_k,t_{k+1}]
  k = MAX(k, p); // usually k>=p, but to avoid round-off
  k = MIN(k, n - 1); // usually k<n, but to avoid round-off
  //////////////////////////
  ASSERT(p <= k && k < n); // because the parameter tnow corresponds to a certain point on the boundary
  //////////////////////////

  const double d = (tnow - t[k]) / dt;

  double u = 0.0; // initialise

  if (p == 2) {

    u += xvec[(k - 0) % (n - p)] * N20(d);
    u += xvec[(k - 1) % (n - p)] * N21(d + 1.0);
    u += xvec[(k - 2) % (n - p)] * N22(d + 2.0);

  } else if (p == 3) {

    u += xvec[(k - 0) % (n - p)] * N30(d);
    u += xvec[(k - 1) % (n - p)] * N31(d + 1.0);
    u += xvec[(k - 2) % (n - p)] * N32(d + 2.0);
    u += xvec[(k - 3) % (n - p)] * N33(d + 3.0);

  } else if (p == 4) {

    u += xvec[(k - 0) % (n - p)] * N40(d);
    u += xvec[(k - 1) % (n - p)] * N41(d + 1.0);
    u += xvec[(k - 2) % (n - p)] * N42(d + 2.0);
    u += xvec[(k - 3) % (n - p)] * N43(d + 3.0);
    u += xvec[(k - 4) % (n - p)] * N44(d + 4.0);

  } else {
    INFO("p=%d is not implemented yet. Exit.\n", p);
    exit(EXIT_FAILURE);
  }

  return u;
}


double2 closedcurve_position(const int p, const int n, const double *t, const double2 *cp, const double tnow)
{
  const double dt = 1.0 / (n + p); // interval for uniform knots
  int k = (int)(tnow / dt); // current point is in the interval [t_k,t_{k+1}]
  k = MAX(k, p); // usually k>=p, but to avoid round-off
  k = MIN(k, n - 1); // usually k<n, but to avoid round-off
  //////////////////////////
  ASSERT(p <= k && k < n); // because the parameter tnow corresponds to a certain point on the boundary
  //////////////////////////

  const double d = (tnow - t[k]) / dt;

  double2 y;
  y.x = y.y = 0.0; // initialise

  if (p == 2) {

    y = add2(y, scale2(N20(d),       cp[k - 0]));
    y = add2(y, scale2(N21(d + 1.0), cp[k - 1]));
    y = add2(y, scale2(N22(d + 2.0), cp[k - 2]));

  } else if (p == 3) {

    y = add2(y, scale2(N30(d),       cp[k - 0]));
    y = add2(y, scale2(N31(d + 1.0), cp[k - 1]));
    y = add2(y, scale2(N32(d + 2.0), cp[k - 2]));
    y = add2(y, scale2(N33(d + 3.0), cp[k - 3]));

  } else if (p == 4) {

    y = add2(y, scale2(N40(d),       cp[k - 0]));
    y = add2(y, scale2(N41(d + 1.0), cp[k - 1]));
    y = add2(y, scale2(N42(d + 2.0), cp[k - 2]));
    y = add2(y, scale2(N43(d + 3.0), cp[k - 3]));
    y = add2(y, scale2(N44(d + 4.0), cp[k - 4]));

  } else {
    INFO("p=%d is not implemented yet. Exit.\n", p);
    exit(EXIT_FAILURE);
  }

  return y;
}


double2 closedcurve_position_derivative(const int p, const int n, const double *t, const double2 *cp, const double tnow)
{
  const double dt = 1.0 / (n + p); // interval for uniform knots
  int k = (int)(tnow / dt); // current point is in the interval [t_k,t_{k+1}]
  k = MAX(k, p); // usually k>=p, but to avoid round-off
  k = MIN(k, n - 1); // usually k<n, but to avoid round-off
  //////////////////////////
  ASSERT(p <= k && k < n); // because the parameter tnow corresponds to a certain point on the boundary
  //////////////////////////

  const double d = (tnow - t[k]) / dt;

  double2 dy;
  dy.x = dy.y = 0.0; // initialise

  if (p == 2) {

    dy = add2(dy, scale2(DN20(d),       cp[k - 0]));
    dy = add2(dy, scale2(DN21(d + 1.0), cp[k - 1]));
    dy = add2(dy, scale2(DN22(d + 2.0), cp[k - 2]));

  } else if (p == 3) {

    dy = add2(dy, scale2(DN30(d),       cp[k - 0]));
    dy = add2(dy, scale2(DN31(d + 1.0), cp[k - 1]));
    dy = add2(dy, scale2(DN32(d + 2.0), cp[k - 2]));
    dy = add2(dy, scale2(DN33(d + 3.0), cp[k - 3]));

  } else if (p == 4) {

    dy = add2(dy, scale2(DN40(d),       cp[k - 0]));
    dy = add2(dy, scale2(DN41(d + 1.0), cp[k - 1]));
    dy = add2(dy, scale2(DN42(d + 2.0), cp[k - 2]));
    dy = add2(dy, scale2(DN43(d + 3.0), cp[k - 3]));
    dy = add2(dy, scale2(DN44(d + 4.0), cp[k - 4]));

  } else {
    INFO("p=%d is not implemented yet. Exit.\n", p);
    exit(EXIT_FAILURE);
  }

  dy = scale2(1.0 / dt, dy); // due to the differentialtion of (tnow-t[j])/dt

  return dy;
}


double2 closedcurve_position_derivative2(const int p, const int n, const double *t, const double2 *cp, const double tnow)
{
  const double dt = 1.0 / (n + p); // interval for uniform knots
  int k = (int)(tnow / dt); // current point is in the interval [t_k,t_{k+1}]
  k = MAX(k, p); // usually k>=p, but to avoid round-off
  k = MIN(k, n - 1); // usually k<n, but to avoid round-off
  //////////////////////////
  ASSERT(p <= k && k < n); // because the parameter tnow corresponds to a certain point on the boundary
  //////////////////////////

  const double d = (tnow - t[k]) / dt;

  double2 ddy;
  ddy.x = ddy.y = 0.0; // initialise

  if (p == 2) {

    ddy = add2(ddy, scale2(DDN20(d),       cp[k - 0]));
    ddy = add2(ddy, scale2(DDN21(d + 1.0), cp[k - 1]));
    ddy = add2(ddy, scale2(DDN22(d + 2.0), cp[k - 2]));

  } else if (p == 3) {

    ddy = add2(ddy, scale2(DDN30(d),       cp[k - 0]));
    ddy = add2(ddy, scale2(DDN31(d + 1.0), cp[k - 1]));
    ddy = add2(ddy, scale2(DDN32(d + 2.0), cp[k - 2]));
    ddy = add2(ddy, scale2(DDN33(d + 3.0), cp[k - 3]));

  } else if (p == 4) {

    ddy = add2(ddy, scale2(DDN40(d),       cp[k - 0]));
    ddy = add2(ddy, scale2(DDN41(d + 1.0), cp[k - 1]));
    ddy = add2(ddy, scale2(DDN42(d + 2.0), cp[k - 2]));
    ddy = add2(ddy, scale2(DDN43(d + 3.0), cp[k - 3]));
    ddy = add2(ddy, scale2(DDN44(d + 4.0), cp[k - 4]));

  } else {
    INFO("p=%d is not implemented yet. Exit.\n", p);
    exit(EXIT_FAILURE);
  }

  ddy = scale2(1.0 / (dt * dt), ddy); // due to the differentialtion of (tnow-t[j])/dt

  return ddy;
}


double2 closedcurve_position_derivative3(const int p, const int n, const double *t, const double2 *cp, const double tnow)
{
  const double dt = 1.0 / (n + p); // interval for uniform knots
  int k = (int)(tnow / dt); // current point is in the interval [t_k,t_{k+1}]
  k = MAX(k, p); // usually k>=p, but to avoid round-off
  k = MIN(k, n - 1); // usually k<n, but to avoid round-off
  //////////////////////////
  ASSERT(p <= k && k < n); // because the parameter tnow corresponds to a certain point on the boundary
  //////////////////////////

  const double d = (tnow - t[k]) / dt;

  double2 dddy;
  dddy.x = dddy.y = 0.0; // initialise

  if (p == 2) {

#if defined(SLOW)
    dddy = add2(dddy, scale2(DDDN20(d),       cp[k - 0]));
    dddy = add2(dddy, scale2(DDDN21(d + 1.0), cp[k - 1]));
    dddy = add2(dddy, scale2(DDDN22(d + 2.0), cp[k - 2]));
#else
    dddy.x = dddy.y = 0.0; // zero
#endif

  } else if (p == 3) {

    dddy = add2(dddy, scale2(DDDN30(d),       cp[k - 0]));
    dddy = add2(dddy, scale2(DDDN31(d + 1.0), cp[k - 1]));
    dddy = add2(dddy, scale2(DDDN32(d + 2.0), cp[k - 2]));
    dddy = add2(dddy, scale2(DDDN33(d + 3.0), cp[k - 3]));

  } else if (p == 4) {

    dddy = add2(dddy, scale2(DDDN40(d),       cp[k - 0]));
    dddy = add2(dddy, scale2(DDDN41(d + 1.0), cp[k - 1]));
    dddy = add2(dddy, scale2(DDDN42(d + 2.0), cp[k - 2]));
    dddy = add2(dddy, scale2(DDDN43(d + 3.0), cp[k - 3]));
    dddy = add2(dddy, scale2(DDDN44(d + 4.0), cp[k - 4]));

  } else {
    INFO("p=%d is not implemented yet. Exit.\n", p);
    exit(EXIT_FAILURE);
  }

  dddy = scale2(1.0 / (dt * dt * dt), dddy); // due to the differentialtion of (tnow-t[j])/dt

  return dddy;
}
