#include "nbie.h"

#if defined(NEWTONCOTES)
#include "newtoncotes.h"
#endif

static double qinc(const int ninc, const double2 x)
{
  double q;

  switch (ninc) {
  case 100:
    abort();
    break;
    
  case 200:
    abort();
    break;
    
  case 300: //This is for test purpose (external problem)
    q = - x.y; // n_2=-x_2 
    break;
    
  default:
    MESG("Unregistered ninc was used. Exit.\n");
    exit(EXIT_FAILURE);
  }

  return q;
}


void mkqinc(const int ncol, const int ninc, const double2 *xc, double *bvec)
{

  if (ninc != NULL_UINC) {

    for (int ix = 0; ix < ncol; ix ++) {

      const double q = qinc(ninc, xc[ix]);

      bvec[ix] += q;

    }

  }

}

void isolap2d_mkqinc(const Isolap2d *p, const int ncol, double *bvec)
{
  const Params *pa = p->params;
  mkqinc(ncol, pa->ninc, p->xc, bvec);
}

/* Numerical integral(p=2) */
double num_integral22(const GaussOne *g, const double a, const double b, const double c, const double d, const double e, const double f, const double x1, const double x2)
{
  double In = 0.0;
  
  const int ng = g->n;
  const double *gx = g->x;
  const double *gw = g->w;
  
  const double Jg = (x2 - x1) / 2.0; // Jacobian

  for (int ig = 0; ig < ng; ig ++) {	    
    const double xnow = ((1.0 - gx[ig]) * x1 + (1.0 + gx[ig]) * x2) / 2.0;
      
    const double kern = ((d * xnow + e) * xnow + f) / ((a * xnow + b) * xnow + c);
    
    In += gw[ig] * Jg * kern;
    
  } // ig
  
  return In;

 }

/* Numerical integral(p=2) */
double num_integral24(const GaussOne *g, const double a, const double b, const double c, const double d, const double e, const double f, const double x1, const double x2)
{
  double In = 0.0;
  
  const int ng = g->n;
  const double *gx = g->x;
  const double *gw = g->w;
  
  const double Jg = (x2 - x1) / 2.0; // Jacobian

  for (int ig = 0; ig < ng; ig ++) {	    
    const double xnow = ((1.0 - gx[ig]) * x1 + (1.0 + gx[ig]) * x2) / 2.0;
      
    const double kern = ((d * xnow + e) * xnow + f) / (((a * xnow + b) * xnow + c) * ((a * xnow + b) * xnow + c));
    
    In += gw[ig] * Jg * kern;
    
  } // ig
  
  return In;

 }

//140318/* Numerical integral(p=3) */
//140318double cbie_integral24(const GaussOne *g, const double2 va, const double2 vb, const double2 vc, const double a, const double b, const double c, const double d, const double e, const double f, const double U0, const double U1, const double U2, const double U3, const double x1, const double x2)
//140318{
//140318  double In = 0.0;
//140318  
//140318  const int ng = g->n;
//140318  const double *gx = g->x;
//140318  const double *gw = g->w;
//140318  
//140318  const double Jg = (x2 - x1) / 2.0; // Jacobian
//140318
//140318  for (int ig = 0; ig < ng; ig ++) {	    
//140318    const double xnow = ((1.0 - gx[ig]) * x1 + (1.0 + gx[ig]) * x2) / 2.0;
//140318    const double rr =  (((a * xnow + b) * xnow + (c + d)) * xnow + e) * xnow + f;
//140318          
//140318    //140225    const double kern = ((va.y * vb.x - va.x * vb.y) / 2.0 + xnow * (va.x * vc.y - va.y * vc.x) / 3.0 + xnow * xnow * (vb.y * vc.x - vb.x * vc.y) / 6.0) / rr * (((U3 * xnow + U2) * xnow + U1) * xnow + U0);
//140318    const double kern = ((va.y * vb.x - va.x * vb.y) / 2.0 - 2.0 * xnow * (va.x * vc.y - va.y * vc.x) / 3.0 + xnow * xnow * (vb.y * vc.x - vb.x * vc.y) / 6.0) / rr * (((U3 * xnow + U2) * xnow + U1) * xnow + U0);
//140318    
//140318    In += gw[ig] * Jg * kern;
//140318
//140318  } // ig
//140318  
//140318  return In;
//140318
//140318 }
//140318
//140318/* Numerical integral(f1,p=3) */
//140318double num_integral54(const GaussOne *g, const double a, const double b, const double c, const double d, const double e, const double f, const double U0, const double U1, const double U2, const double U3, const double x1, const double x2)
//140318{
//140318  double In = 0.0;
//140318  
//140318  const int ng = g->n;
//140318  const double *gx = g->x;
//140318  const double *gw = g->w;
//140318  
//140318  const double Jg = (x2 - x1) / 2.0; // Jacobian
//140318
//140318  for (int ig = 0; ig < ng; ig ++) {	    
//140318    const double xnow = ((1.0 - gx[ig]) * x1 + (1.0 + gx[ig]) * x2) / 2.0;
//140318    const double rr =  (((a * xnow + b) * xnow + (c + d)) * xnow + e) * xnow + f;
//140318
//140318    const double kern = ((a * xnow + b) * xnow + c - d / 2.0) / rr * (((U3 * xnow + U2) * xnow + U1) * xnow + U0);
//140318    
//140318    In += gw[ig] * Jg * kern;
//140318    
//140318  } // ig
//140318  
//140318  return In;
//140318
//140318 }
//140318
//140318/* Numerical integral(f2,p=3) */
//140318double num_integral68(const GaussOne *g, const double2 va, const double2 vb, const double2 vc, const double a, const double b, const double c, const double d, const double e, const double f, const double U0, const double U1, const double U2, const double U3, const double x1, const double x2)
//140318{
//140318  double In = 0.0;
//140318  
//140318  const int ng = g->n;
//140318  const double *gx = g->x;
//140318  const double *gw = g->w;
//140318  
//140318  const double Jg = (x2 - x1) / 2.0; // Jacobian
//140318
//140318  for (int ig = 0; ig < ng; ig ++) {	    
//140318    const double xnow = ((1.0 - gx[ig]) * x1 + (1.0 + gx[ig]) * x2) / 2.0;
//140318    const double rr =  (((a * xnow + b) * xnow + (c + d)) * xnow + e) * xnow + f;
//140318      
//140318    //140225    const double kern = ((va.x * vb.y - va.y * vb.x) / 2.0 + xnow * (va.x * vc.y - va.y * vc.x) / 3.0) * ((va.y * vb.x - va.x * vb.y) / 2.0 + xnow * (va.x * vc.y - va.y * vc.x) / 3.0 + xnow * xnow * (vb.y * vc.x - vb.x * vc.y) / 6.0) / (rr * rr) * (((U3 * xnow + U2) * xnow + U1) * xnow + U0);
//140318    const double kern = ((va.y * vb.x - va.x * vb.y) / 2.0 - 2.0 * xnow * (va.x * vc.y - va.y * vc.x) / 3.0 + xnow * xnow * (vb.y * vc.x - vb.x * vc.y) / 6.0) * ((va.x * vb.y - va.y * vb.x) / 2.0 + xnow * (va.x * vc.y - va.y * vc.x) / 3.0) / (rr * rr) * (((U3 * xnow + U2) * xnow + U1) * xnow + U0);
//140318    
//140318    In += gw[ig] * Jg * kern;
//140318    
//140318  } // ig
//140318  
//140318  return In;
//140318
//140318 }



/* integral of f1 or CBIE(p=2) */
double integrate1(const double a, const double b, const double c, const double d, const double e, const double f, const double x)
{
  ASSERT(a >= 0.0); // if no round-off error
  ASSERT(c > 0.0); // if no round-off error

  const double D = b * b - 4.0 * a * c; // discriminant
  ASSERT(D <= 0.0); // if no round-off error

  double In;
  
  if (a == 0.0) {
    //  if (SQUARE(a) <= SQUARE(SMALL_NUMBER)) {

    if (b == 0.0) { // a=0 & b=0
    //    if (SQUARE(b) <= SQUARE(SMALL_NUMBER)) { // a=0 & b=0

      In = x * (f + x * (e / 2.0 + x * d / 3.0)) / c;

    } else { // a=0 & b!=0

      In = (b * x * (- 2.0 * c * d + 2.0 * b * e + b * d * x) + 2.0 * (c * c * d - b * c * e + b * b * f) * log(fabs(c + b * x))) / (2.0 * b * b * b); 

    }
    
  } else {    
    
    if (D == 0.0) {
      //    if (SQUARE(D) <= SQUARE(SMALL_NUMBER)) {
			       
      //140225      In = 1.0 / (a * a * (b + 2.0 * a * x)) * (a * (b * e - 2.0 * a * f + 2.0 * b * d * x + 2.0 * a * d * x * x) - (b * d - a * e) * (b + 2.0 * a * x) * log(fabs(b + 2.0 * a * x)));
      
      In = 1.0 / (a * a) * (1.0 / (b + 2.0 * a * x) * (a * (b * e - 2.0 * a * f + 2.0 * b * d * x + 2.0 * a * d * x * x)) - (b * d - a * e) * log(fabs(b + 2.0 * a * x)));

      } else {
	
	In = 1.0 / (2.0 * a * a) * (2.0 * a * d * x + (2.0 * (b * b * d - a * b * e + 2.0 * a * (- c * d + a * f)) * atan((b + 2.0 * a * x) / sqrt(- D))) / sqrt(- D) + (- b * d + a * e) * log(fabs(c + x * (b + a * x))));
	
    }                
            
  }

  return In;
  
}

/* integral of f2(p=2) */
double integrate2(const double a, const double b, const double c, const double d, const double e, const double f, const double x)
{
  ASSERT(a >= 0.0); // if no round-off error
  ASSERT(c > 0.0); // if no round-off error

  const double D = b * b - 4.0 * a * c; // discriminant
  ASSERT(D <= 0.0); // if no round-off error

  double In;
  
  if (a == 0) {
    
    In = (b * d * x - (c * c * d - b * c * e + b * b * f) / (c + b * x) + (- 2.0 * c * d + b * e) * log(fabs(c + b * x))) / (b * b * b);
      
  }
  
  else {
  
    if (D == 0) {
			       
      In = - (2.0 * (b * b * d + a * b * (e + 6.0 * d * x) + 2.0 * a * a * (2.0 * f + 3.0 * x * (e + 2.0 * d * x)))) / (3.0 * a * (b + 2.0 * a * x) * (b + 2.0 * a * x) * (b + 2.0 * a * x));  
	      
    } else {
	
      In = (b * b * d * x + 2.0 * a * (a * f * x - c * (e + d * x)) + b * (c * d + a * (f - e * x))) / (a * (- D) * (c + x * (b + a * x))) - (2.0 * (- 2.0 * c * d + b * e - 2.0 * a * f) * atan((b + 2.0 * a * x) / (sqrt(- D)))) / (- D * sqrt(- D)); 
	
    }
      
  }
  
  return In;

}



static double sekibun1plus(const double a, const double b, const double c, const double d, const double e, const double f)
{
  /*
    \int_0^1 \frac{d*x^2+e*x+f}{a*x^2+b*x+c} dx
    where a=dot2(vb,vb)/4, b=dot2(va,vb), and c=dot2(va,va)
  */

  ASSERT(a >= 0.0); // if no round-off error
  ASSERT(c > 0.0); // if no round-off error

  double In;

  //  if (SQUARE(a) <= SQUARE(SMALL_NUMBER)) { // if a=0, then b=0
  if (a == 0.0) { // if a=0, then b=0

    In = (d / 3.0 + e / 2.0 + f) / c;

  } else {
    
    //    if (SQUARE(b) <= SQUARE(SMALL_NUMBER)) { // a>0, b=0
    if (b == 0.0) { // a>0, b=0

      //      In = 1.0 / a * (d + (- c * d + a * f) / sqrt(a * c) * atan(sqrt(a / c)) + e * atanh(a / (a + 2.0 * c)));
      In = 1.0 / a * (d + (- c * d + a * f) / sqrt(a * c) * atan(sqrt(a / c)) + e * log((a + c) / c) / 2.0); // atanh(x)=log((1+x)/(1-x)) if -1<x<1

    } else {

      const double D = b * b - 4.0 * a * c; // discriminant
      ASSERT(D <= 0.0); // if no round-off error

      //      if (SQUARE(D) <= SQUARE(SMALL_NUMBER)) { // a>0, b!=0, D=0
      if (D == 0.0) { // a>0, b!=0, D=0

	if (b > 0 || b + 2.0 * a < 0) {

	  In = 2.0 / (a * a * b * (2.0 * a + b)) * (a * b * ((a + b) * d - a * e) + 2.0 * a * a * a * f + b * (2.0 * a + b) * (- b * d + a * e) * atanh(a / (a + b)));
	  //	  In = 2.0 / (a * a * b * (2.0 * a + b)) * (a * b * ((a + b) * d - a * e) + 2.0 * a * a * a * f + b * (2.0 * a + b) * (- b * d + a * e) * log((2.0 * a + b) / b) / 2.0);

	} else { // must be impossible;

	  abort();

	}

      } else { // a>0, b!=0, D<0

	In = 1.0 / (a * a) * ((b * b * d - a * b * e + 2.0 * a * (- c * d + a * f)) / sqrt(- D) * (atan((2.0 * a + b) / sqrt(- D)) - atan(b / sqrt(- D))) + a * d + (b * d - a * e) / 2.0 * log(c / (a + b + c)));

      }

    }

  }

  return In;
  
}


static double sekibun1minus(const double a, const double b, const double c, const double d, const double e, const double f)
{
  /*
    \int_{-1}^0 \frac{d*x^2+e*x+f}{a*x^2+b*x+c} dx
    where a=dot2(vb,vb)/4, b=dot2(va,vb), and c=dot2(va,va)
  */

  ASSERT(a >= 0.0); // if no round-off error
  ASSERT(c > 0.0); // if no round-off error

  double In;

  //  if (SQUARE(a) <= SQUARE(SMALL_NUMBER)) { // if a=0, then b=0
  if (a == 0.0) { // if a=0, then b=0

    In = (d / 3.0 - e / 2.0 + f) / c;

  } else {
    
    //    if (SQUARE(b) <= SQUARE(SMALL_NUMBER)) { // a>0, b=0
    if (b == 0.0) { // a>0, b=0

      In = 1.0 / a * (d + (- c * d + a * f) / sqrt(a * c) * atan(sqrt(a / c)) + e * log(c / (a + c)) / 2.0); // atanh(x)=log((1+x)/(1-x)) if -1<x<1

    } else {

      const double D = b * b - 4.0 * a * c; // discriminant
      ASSERT(D <= 0.0); // if no round-off error

      //      if (SQUARE(D) <= SQUARE(SMALL_NUMBER)) { // a>0, b!=0, D=0
      if (D == 0.0) { // a>0, b!=0, D=0

	if (b < 0 || b - 2.0 * a > 0) {

	  In = 2.0 / (a * a * b * (2.0 * a - b)) * (a * b * ((a - b) * d + a * e) - 2.0 * a * a * a * f - b * (- 2.0 * a + b) * (b * d - a * e) * atanh(a / (a - b)));

	} else { // must be impossible;

	  abort();

	}

      } else { // a>0, b!=0, D<0

	In = 1.0 / (a * a) * ((b * b * d - a * b * e + 2.0 * a * (- c * d + a * f)) / sqrt(- D) * (atan(b / sqrt(- D)) - atan((- 2.0 * a + b) / sqrt(- D))) + a * d + (a * e - b * d) / 2.0 * log(c / (a - b + c)));

      }

    }

  }

  return In;
  
}



static double num_sekibun1(const GaussOne *g, const double a, const double b, const double c, const double d, const double e, const double f, const double x1, const double x2)
{
  double In = 0.0;
  
  const int ng = g->n;
  const double *gx = g->x;
  const double *gw = g->w;
  
  const double Jg = (x2 - x1) / 2.0; // Jacobian
  
  for (int ig = 0; ig < ng; ig ++) {	    
    const double xnow = ((1.0 - gx[ig]) * x1 + (1.0 + gx[ig]) * x2) / 2.0;
    
    const double kern = ((d * xnow + e) * xnow + f) / ((a * xnow + b) * xnow + c);
    
    In += gw[ig] * Jg * kern;
    
  } // ig
  
  return In;

}


#ifdef MYINTDE

#include "myintde2.h"

double myfunc(double tnow, Intde *intde)
{
  const int p = intde->params->p;
  const int n = intde->params->n;
  const double *t = intde->params->t;
  const double2 *cp = intde->params->cp;
  const double2 x = intde->params->x;
  const double2 nx = intde->params->nx;
  const double *xvec = intde->params->xvec;
  const int joff = intde->params->joff;
  const int k = intde->params->k;
  
  const double2 y = closedcurve_position(p, n, t, cp, tnow);
  const double2 dy = closedcurve_position_derivative(p, n, t, cp, tnow);
  const double2 vr = sub2(x, y);
  const double rr = dot2(vr, vr);
  
#if defined(CBIE)
  const double kern = (vr.x * dy.y - vr.y * dy.x) / rr;
#else
  const double kern = ((nx.x * dy.y - nx.y * dy.x) - 2.0 * (vr.x * nx.x + vr.y * nx.y) * (vr.x * dy.y - vr.y * dy.x) / rr) / rr;
#endif
  
  double bsum = 0.0;
  for (int l = 0; l <= p; l ++) { // then, k-p<=k-l<=k, i.e., j<=k-l<=j+p
    const double bsp = closedcurve_basis(p, n, t, k - l, tnow);
    bsum += bsp * xvec[joff + (k - l) % (n - p)];
  }
  
  return kern * bsum;
}

#endif




static double sekibun_p2(const GaussOne *g, const double2 va, const double2 vb, const double U0, const double U1, const double U2, const double x1, const double x2)
{

  const double a = dot2(vb, vb) / 4.0;
  const double b = dot2(va, vb);
  const double c = dot2(va, va);

  double In;

#if defined(CBIE)

  //140318  const double cf0 = - (va.x * vb.y - va.y * vb.x) / 2.0;
  const double cf0 = - vdot2(va, vb) / 2.0;
#if defined(NUM)
  In = cf0 * num_sekibun1(g, a, b, c, U2, U1, U0, x1, x2); // N21, p10
#else
  //140318#if(1)
  In = cf0 * (integrate1(a, b, c, U2, U1, U0, x2) - integrate1(a, b, c, U2, U1, U0, x1));
  //140318#else
  //140318  In = cf0 * sekibun1plus(a, b, c, U2, U1, U0); // N21, p10; same as the above
  //140318#endif
#endif

#else

  const double cf1 = - a / norm2(va); // common factor of f1
  //140318  const double cf2 = - (va.x * vb.y - va.y * vb.x) * (va.x * vb.y - va.y * vb.x) / (2.0 * norm2(va)); // common factor of f2
  const double cf2 = - vdot2(va, vb) * (va.x * vb.y - va.y * vb.x) / (2.0 * norm2(va)); // common factor of f2
#if defined(NUM)
  In = cf1 * num_integral22(g, a, b, c, U2, U1, U0, x1, x2) - cf2 * num_integral24(g, a, b, c, U2, U1, U0, x1, x2);
#else
  In = cf1 * (integrate1(a, b, c, U2, U1, U0, x2) - integrate1(a, b, c, U2, U1, U0, x1)) - cf2 * (integrate2(a, b, c, U2, U1, U0, x2) - integrate2(a, b, c, U2, U1, U0, x1));
#endif

  In += 1.0 / norm2(va) * ((- 1.0) * U0 + 0.0 * U1 + (x2 - x1) * U2 / 1.0);

#endif

  return In;

 }


static double sekibun_p3(const GaussOne *g, const double2 va, const double2 vb, const double2 vc, const double U0, const double U1, const double U2, const double U3, const double x1, const double x2)
{
  double In = 0.0;
  
  const int ng = g->n;
  const double *gx = g->x;
  const double *gw = g->w;
  
  const double Jg = (x2 - x1) / 2.0; // Jacobian

  for (int ig = 0; ig < ng; ig ++) {	    
    const double xnow = ((1.0 - gx[ig]) * x1 + (1.0 + gx[ig]) * x2) / 2.0;
    const double2 vp = add2(scale2(1.0 / 2.0, vb), scale2(xnow / 3.0, vc)); // b/2+c/3*x
    const double2 vq = add2(vb, scale2(xnow, vc)); // b+c*x
    const double2 vr = add2(va, scale2(xnow, vp)); // -r/x:=a+b/2*x+c/3*x^2
    const double rr = dot2(vr, vr); // r^2/x^2
#if defined(CBIE)
    const double kern = (- vdot2(va, vq) - vdot2(vp, va) - xnow * vdot2(vp, vq)) / rr;
#else
    const double kern = ((dot2(va, scale2(1.0 / 3.0, vc)) - dot2(vp, vp)) - 2.0 * vdot2(vp, va) * (vdot2(va, vq) + vdot2(vp, va) + vdot2(vp, vq)) / rr) / rr;
#endif
    const double u = U0 + xnow * (U1 + xnow * (U2 + xnow * U3));
    In += gw[ig] * Jg * kern * u;
  }
  
#if defined(CBIE)
#else
  In += (- 1.0) * U0 + 0.0 * U1 + (x2 - x1) * (U2 / 1.0 + (x2 - x1) * U3 / 2.0);
  In /= norm2(va);
#endif

  return In;

 }


static double sekibun_p4(const GaussOne *g, const double2 va, const double2 vb, const double2 vc, const double2 vd, const double U0, const double U1, const double U2, const double U3, const double U4, const double x1, const double x2)
{
  double In = 0.0;

  const int ng = g->n;
  const double *gx = g->x;
  const double *gw = g->w;
  
  const double Jg = (x2 - x1) / 2.0; // Jacobian

  for (int ig = 0; ig < ng; ig ++) {	    
    const double xnow = ((1.0 - gx[ig]) * x1 + (1.0 + gx[ig]) * x2) / 2.0;
    const double2 vp = add2(scale2(1.0 / 2.0, vb), add2(scale2(xnow / 3.0, vc), scale2(xnow * xnow / 4.0, vd))); // b/2+c/3*x+d/4*x^2
    const double2 vq = add2(vb, add2(scale2(xnow, vc), scale2(xnow * xnow, vd))); // b+c*x+d*x^2
    //    const double2 vr = scale2(- 1.0, add2(va, scale2(xnow, vp))); // r/x:=-(a+b/2*x+c/3*x^2+d/4*x^3)
    const double2 vr = add2(va, scale2(xnow, vp)); // -r/x:=a+b/2*x+c/3*x^2+d/4*x^3
    const double rr = dot2(vr, vr); // r^2/x^2
#if defined(CBIE)
    //    const double kern = ((va.y * vq.x - va.x * vq.y) + (vp.y * va.x - vp.x * va.y) + xnow * (vp.y * vq.x - vp.x * vq.y)) / rr; // (r.n/x^2)/(r^2/x^2)
    const double kern = (- vdot2(va, vq) - vdot2(vp, va) - xnow * vdot2(vp, vq)) / rr; // (r.n/x^2)/(r^2/x^2)
#else
    //    double kern = - coef * (dot2(vp, vp) - dot2(va, add2(scale2(1.0 / 3.0, vc), scale2(xnow / 2.0, vd)))) / rr; // f1
    //    kern -= coef * 2.0 * vdot2(vp, va) * (vdot2(va, vq) + vdot2(vp, va) + vdot2(vp, vq)) / (rr * rr); // f2
    const double kern = ((dot2(va, add2(scale2(1.0 / 3.0, vc), scale2(xnow / 2.0, vd))) - dot2(vp, vp)) - 2.0 * vdot2(vp, va) * (vdot2(va, vq) + vdot2(vp, va) + vdot2(vp, vq)) / rr) / rr;
#endif
    const double u = U0 + xnow * (U1 + xnow * (U2 + xnow * (U3 + xnow * U4)));
    In += gw[ig] * Jg * kern * u;
  }
  
#if defined(CBIE)
#else
  //  In += U0 * (- 1.0) + U1 * 0.0 + (x2 - x1) * (U2 / 1.0 + (x2 - x1) * (U3 / 2.0 + (x2 - x1) * U4 / 3.0));
  In += (- 1.0) * U0 + 0.0 * U1 + (x2 - x1) * (U2 / 1.0 + (x2 - x1) * (U3 / 2.0 + (x2 - x1) * U4 / 3.0));
  In /= norm2(va);
#endif

  return In;

 }




//140325void isolap2d_direct_nbie(const GaussOne *g, const int p, const int n, const double *t, const double2 *cp, const double *tc, const double2 *xc, const int icc, const int jcc, double *dval, const double *xvec, const int ioff, const int joff, const int i, const int pi, const int ni, const double *ti, const double2 *cpi)
void isolap2d_direct_nbie(const GaussOne *g, const int p, const int n, const double *t, const double2 *cp, const int icc, const int jcc, double *dval, const double *xvec, const int ioff, const int joff, const int i, const double2 x, const double2 nx, const int j)
{  
  //140325  /* Global index of collocation point */
  //140325  const int ix = ioff + i; // i is local
  //140325
  //140325  /* Collocation point x and its derivative dx/dt */
  //140325  const double2 x = xc[ix];
  //140325  //140225  const double2 dx = closedcurve_position_derivative(p, n, t, cp, tc[ix]);
  //140325  const double2 dx = closedcurve_position_derivative(pi, ni, ti, cpi, tc[ix]);
  //140325
  //140325  /* Compute the outward normal vector at the collocation point */
  //140325  double2 nx;
  //140325  nx.x = dx.y / norm2(dx); // nx1
  //140325  nx.y = - dx.x / norm2(dx); // nx2
  
  /* Set up for Guassian quadrature */
  const int ng = g->n;
  const double *gx = g->x;
  const double *gw = g->w;
  
#ifdef MYINTDE
  /* Set up for DE formula */
  Intde *intde = myintde_malloc(); // this can be moved to matvec_nbie
  myintde_init(intde); // this can be moved to matvec_nbie
  IntdeParams intdeparams;
  intde->params = &intdeparams;
  intdeparams.p = p;
  intdeparams.n = n;
  intdeparams.t = t;
  intdeparams.cp = cp;
  intdeparams.x = x;
  intdeparams.nx = nx;
  intdeparams.xvec = xvec;
  intdeparams.joff = joff;
#endif

  //140325  /* Initialise */
  //140325  *dval = 0.0;
  //140325
  //140325  /* Loop over sections; I_j=[t_{j+p},t_{j+p+1}) or [t_k,t_{k+1}), where 0<=j<n-p */
  //140325  for (int j = 0; j < n - p; j ++) {

  const int k = j + p; // p<=k<n

  /* Integral over I_j for x_i:=y(t_{i+p}) */
  if (i == j && icc == jcc) { // Singular at the lefthand-side end of I_j or I_i
    
    const double x1 = 0.0; // lower bound of X
    const double x2 = 1.0; // upper bound of X
    
    if (p == 2) {
      
      /* Compute vectors a and b (where k-1>=1 and k-2>=0) */
      const double2 va = sub2(cp[k - 1], cp[k - 2]);
      const double2 vb = sub2(add2(cp[k], cp[k - 2]), scale2(2.0, cp[k - 1]));
      
      const double u0 = xvec[joff + (k - 0) % (n - p)];
      const double u1 = xvec[joff + (k - 1) % (n - p)];
      const double u2 = xvec[joff + (k - 2) % (n - p)];
      
      const double U0 = (u1 + u2) / 2.0;
      const double U1 = u1 - u2;
      const double U2 = (u0 - 2.0 * u1 + u2) / 2.0;
      
      //140325      *dval += sekibun_p2(g, va, vb, U0, U1, U2, x1, x2);
      *dval = sekibun_p2(g, va, vb, U0, U1, U2, x1, x2);
	
    } else if (p == 3) {
      
      /* Compute vector a, b, and c (p=3); See Eq.(B.5) */
      const double2 va = scale2(0.5, sub2(cp[k - 1], cp[k - 3]));
      const double2 vb = sub2(add2(cp[k - 1], cp[k - 3]), scale2(2.0, cp[k - 2]));
      const double2 vc = sub2(scale2(0.5, sub2(cp[k], cp[k - 3])), scale2(1.5, sub2(cp[k - 1], cp[k - 2])));       
      
      /* Compute potential (p=3); See Eq.(B.11) */      
      const double u0 = xvec[joff + (k - 0) % (n - p)];
      const double u1 = xvec[joff + (k - 1) % (n - p)];
      const double u2 = xvec[joff + (k - 2) % (n - p)];
      const double u3 = xvec[joff + (k - 3) % (n - p)];
      
      const double U0 = (u1 + 4.0 * u2 + u3) / 6.0;
      const double U1 = (u1 - u3) / 2.0;
      const double U2 = (u1 - 2.0 * u2 + u3) / 2.0;
      const double U3 = (u0 - 3.0 * u1 + 3.0 * u2 - u3) / 6.0;
      
      //140325      *dval += sekibun_p3(g, va, vb, vc, U0, U1, U2, U3, x1, x2);
      *dval = sekibun_p3(g, va, vb, vc, U0, U1, U2, U3, x1, x2);

    } else if (p == 4) {
      
      const double2 va = scale2(1.0 / 24.0, add2(scale2(4.0, cp[k - 1]), add2(scale2(12.0, cp[k - 2]), add2(scale2(- 12.0, cp[k - 3]), scale2(- 4.0, cp[k - 4])))));
      const double2 vb = scale2(1.0 / 24.0, add2(scale2(12.0, cp[k - 1]), add2(scale2(- 12.0, cp[k - 2]), add2(scale2(- 12.0, cp[k - 3]), scale2(12.0, cp[k - 4])))));
      const double2 vc = scale2(1.0 / 24.0, add2(scale2(12.0, cp[k - 1]), add2(scale2(- 36.0, cp[k - 2]), add2(scale2(36.0, cp[k - 3]), scale2(- 12.0, cp[k - 4])))));
      const double2 vd = scale2(1.0 / 24.0, add2(scale2(4.0, cp[k - 0]), add2(scale2(- 16.0, cp[k - 1]), add2(scale2(24.0, cp[k - 2]), add2(scale2(- 16.0, cp[k - 3]), scale2(4.0, cp[k - 4]))))));
      
      const double u0 = xvec[joff + (k - 0) % (n - p)];
      const double u1 = xvec[joff + (k - 1) % (n - p)];
      const double u2 = xvec[joff + (k - 2) % (n - p)];
      const double u3 = xvec[joff + (k - 3) % (n - p)];
      const double u4 = xvec[joff + (k - 4) % (n - p)];
      
      const double U0 = (           u1 + 11.0 * u2 + 11.0 * u3 +       u4) / 24.0;
      const double U1 = (     4.0 * u1 + 12.0 * u2 - 12.0 * u3 - 4.0 * u4) / 24.0;
      const double U2 = (     6.0 * u1 -  6.0 * u2 -  6.0 * u3 + 6.0 * u4) / 24.0;
      const double U3 = (     4.0 * u1 - 12.0 * u2 + 12.0 * u3 - 4.0 * u4) / 24.0;
      const double U4 = (u0 - 4.0 * u1 +  6.0 * u2 -  4.0 * u3 +       u4) / 24.0;
      
      //140325      *dval += sekibun_p4(g, va, vb, vc, vd, U0, U1, U2, U3, U4, x1, x2);
      *dval = sekibun_p4(g, va, vb, vc, vd, U0, U1, U2, U3, U4, x1, x2);

    } else { 
      
      /* p>=5 is not available */
      abort();
      
    }
    
    //140225    } else if ((i == j + 1 || (i == 0 && j == n - p - 1)) && icc == jcc) {
  } else if (i == (j + 1) % (n - p) && icc == jcc) { // Singular at the righthand-side of I_j or I_{i-1}=[t_{i+p-1},t_{i+p})
      
    const double x1 = - 1.0; // lower bound of X
    const double x2 = 0.0; // upper bound of X
    
    if (p == 2) {
      
      /* Compute vectors a and b (where k-1>=1 and k-2>=0) */
      const double2 va = sub2(cp[k], cp[k - 1]);
      const double2 vb = sub2(add2(cp[k], cp[k - 2]), scale2(2.0, cp[k - 1]));
      
      /* Compute the coefficients of the enumerator (potential) for the kernel f1 */      
      const double u0 = xvec[joff + (k - 0) % (n - p)];
      const double u1 = xvec[joff + (k - 1) % (n - p)];
      const double u2 = xvec[joff + (k - 2) % (n - p)];
      
      const double U0 = (u0 + u1) / 2.0;
      const double U1 = u0 - u1;
      const double U2 = (u0 - 2.0 * u1 + u2) / 2.0;
      
      //140325      *dval += sekibun_p2(g, va, vb, U0, U1, U2, x1, x2);
      *dval = sekibun_p2(g, va, vb, U0, U1, U2, x1, x2);

    } else if (p == 3) { 
      
      /* Compute vector a, b, and c (p=3); See Eq.(B.23) */
      const double2 va = scale2(0.5, sub2(cp[k], cp[k - 2]));
      const double2 vb = sub2(add2(cp[k], cp[k - 2]), scale2(2.0, cp[k - 1]));
      const double2 vc = sub2(scale2(0.5, sub2(cp[k], cp[k - 3])), scale2(1.5, sub2(cp[k - 1], cp[k - 2])));
      
      /* Compute potential (p=3); See Eq.(B.29) */
      const double u0 = xvec[joff + (k - 0) % (n - p)];
      const double u1 = xvec[joff + (k - 1) % (n - p)];
      const double u2 = xvec[joff + (k - 2) % (n - p)];
      const double u3 = xvec[joff + (k - 3) % (n - p)];
      
      const double U0 = (u0 + 4.0 * u1 + u2) / 6.0;
      const double U1 = (u0 - u2) / 2.0;
      const double U2 = (u0 - 2.0 * u1 + u2) / 2.0;
      const double U3 = (u0 - 3.0 * u1 + 3.0 * u2 - u3) / 6.0;
      
      //140325      *dval += sekibun_p3(g, va, vb, vc, U0, U1, U2, U3, x1, x2);
      *dval = sekibun_p3(g, va, vb, vc, U0, U1, U2, U3, x1, x2);
      
    } else if (p == 4) {
      
      const double2 va = scale2(1.0 / 24.0, add2(scale2(4.0, cp[k - 0]), add2(scale2(12.0, cp[k - 1]), add2(scale2(- 12.0, cp[k - 2]), scale2(- 4.0, cp[k - 3])))));
      const double2 vb = scale2(1.0 / 24.0, add2(scale2(12.0, cp[k - 0]), add2(scale2(- 12.0, cp[k - 1]), add2(scale2(- 12.0, cp[k - 2]), scale2(12.0, cp[k - 3])))));
      const double2 vc = scale2(1.0 / 24.0, add2(scale2(12.0, cp[k - 0]), add2(scale2(- 36.0, cp[k - 1]), add2(scale2(36.0, cp[k - 2]), scale2(- 12.0, cp[k - 3])))));
      const double2 vd = scale2(1.0 / 24.0, add2(scale2(4.0, cp[k - 0]), add2(scale2(- 16.0, cp[k - 1]), add2(scale2(24.0, cp[k - 2]), add2(scale2(- 16.0, cp[k - 3]), scale2(4.0, cp[k - 4]))))));
      
      const double u0 = xvec[joff + (k - 0) % (n - p)];
      const double u1 = xvec[joff + (k - 1) % (n - p)];
      const double u2 = xvec[joff + (k - 2) % (n - p)];
      const double u3 = xvec[joff + (k - 3) % (n - p)];
      const double u4 = xvec[joff + (k - 4) % (n - p)];
      
      const double U0 = (      u0 + 11.0 * u1 + 11.0 * u2 +       u3) / 24.0;
      const double U1 = (4.0 * u0 + 12.0 * u1 - 12.0 * u2 - 4.0 * u3) / 24.0;
      const double U2 = (6.0 * u0 -  6.0 * u1 -  6.0 * u2 + 6.0 * u3) / 24.0;
      const double U3 = (4.0 * u0 - 12.0 * u1 + 12.0 * u2 - 4.0 * u3) / 24.0;
      const double U4 = (      u0 -  4.0 * u1 +  6.0 * u2 - 4.0 * u3 + u4) / 24.0;
      
      //140325      *dval += sekibun_p4(g, va, vb, vc, vd, U0, U1, U2, U3, U4, x1, x2);
      *dval = sekibun_p4(g, va, vb, vc, vd, U0, U1, U2, U3, U4, x1, x2);

    } else { 
      
      /* p>=5 is not available */
      abort();
      
    }
    
    
  } else { // Non-singular integral; See Eq.(2.23) in Sakabe's thesis (which looks to miss the factor '2' in the enumerator of the second term in the RHS)
    
#if defined(NEWTONCOTES)
    
    /* Set up the ng-points Newton-Cotes formula */
    Newtoncotes *nc;
    newtoncotes_load(ng, &nc);
    
    double sum = 0.0; // integral over I_j for x_i
    
    for (int ig = 0; ig < ng; ig ++) {
      const double tnow = (t[k + 1] - t[k]) * nc->h * ig + t[k];
      const double2 y = closedcurve_position(p, n, t, cp, tnow);
      const double2 dy = closedcurve_position_derivative(p, n, t, cp, tnow);
      const double2 vr = sub2(x, y);
      const double rr = dot2(vr, vr);
#if defined(CBIE)
      const double kern = (vr.x * dy.y - vr.y * dy.x) / rr;
#else
      const double kern = ((nx.x * dy.y - nx.y * dy.x) - 2.0 * (vr.x * nx.x + vr.y * nx.y) * (vr.x * dy.y - vr.y * dy.x) / rr) / rr;
#endif
      double bsum = 0.0;
      for (int l = 0; l <= p; l ++) { // then, k-p<=k-l<=k, i.e., j<=k-l<=j+p
	const double bsp = closedcurve_basis(p, n, t, k - l, tnow);
	bsum += bsp * xvec[joff + (k - l) % (n - p)];
      }
      sum += nc->w[ig] * kern * bsum;
    }
    sum *= (t[k + 1] - t[k]) * nc->h;
    
    newtoncotes_free(nc);
    
#elif defined(MYINTDE)
    
    double sum = 0.0; // integral over I_j for x_i
    
    intdeparams.k = k; // k depends on j
    double val, err;
    myintde(intde, &myfunc, t[k], t[k + 1], &val, &err);
    
    sum += val;
    
#elif defined(ADAPTIVE_GAUSS)
    
#if (!defined(ADAPTIVE_GAUSS_EPS) || !defined(ADAPTIVE_GAUSS_NMIN) || !defined(ADAPTIVE_GAUSS_NMAX))
#error Missing.
#endif
    
    const double Jg = (t[k + 1] - t[k]) / 2.0;
    
    int ngtmp = ADAPTIVE_GAUSS_NMIN, converge;
    double sumprev = 0.0;
    
    do {
      
      GaussOne *gtmp;
      GaussOne_load(ngtmp, &gtmp);

      double sumtmp = 0.0;
      for (int ig = 0; ig < ngtmp; ig ++) {
	const double tnow = ((1.0 - gtmp->x[ig]) * t[k] + (1.0 + gtmp->x[ig]) * t[k + 1]) / 2.0;
	const double2 y = closedcurve_position(p, n, t, cp, tnow);
	const double2 dy = closedcurve_position_derivative(p, n, t, cp, tnow);
	const double2 vr = sub2(x, y);
	const double rr = dot2(vr, vr);
#if defined(CBIE)
	const double kern = (vr.x * dy.y - vr.y * dy.x) / rr;
#else
	const double kern = ((nx.x * dy.y - nx.y * dy.x) - 2.0 * (vr.x * nx.x + vr.y * nx.y) * (vr.x * dy.y - vr.y * dy.x) / rr) / rr;
#endif
	double uval = 0.0;
	for (int l = 0; l <= p; l ++) { // then, k-p<=k-l<=k, i.e., j<=k-l<=j+p
	  const double bsp = closedcurve_basis(p, n, t, k - l, tnow);
	  uval += bsp * xvec[joff + (k - l) % (n - p)];
	}
	sumtmp += gtmp->w[ig] * Jg * kern * uval;
      } // ig

      /* Check to go next */
      if (sumprev == 0.0) {
	converge = (sumtmp == 0.0 ? 1 : 0);
      } else {
	converge = (fabs((sumtmp - sumprev) / sumprev) <= ADAPTIVE_GAUSS_EPS ? 1 : 0);
      }

      /* Update */
      sumprev = sumtmp;
      GaussOne_free(gtmp);
      ngtmp ++;

    } while (converge == 0 && ngtmp <= ADAPTIVE_GAUSS_NMAX);

    double sum = sumprev; // integral over I_j for x_i

#else

    double sum = 0.0; // integral over I_j for x_i

    //140325    /* Compute the ends of I_j */
    //140325    const double tminus = t[k];
    //140325    const double tplus = t[k + 1];

    /* Jacobian wrt Gaussian quadrature (this is not Jacobian wrt the curve parameter) */
    //140325    const double Jg = (tplus - tminus) / 2.0; // (t_{k+1}-t_{k})/2
    const double Jg = (t[k + 1] - t[k]) / 2.0;

    for (int ig = 0; ig < ng; ig ++) {
      /* Parameter corresponding to the ig th abscissa */
      //140325      const double tnow = ((1.0 - gx[ig]) * tminus + (1.0 + gx[ig]) * tplus) / 2.0;
      const double tnow = ((1.0 - gx[ig]) * t[k] + (1.0 + gx[ig]) * t[k + 1]) / 2.0;
      /* Integral point y and its derivative dy/dt */	      
      const double2 y = closedcurve_position(p, n, t, cp, tnow);
      const double2 dy = closedcurve_position_derivative(p, n, t, cp, tnow);
      /* Compute r:=|x-y| */
      const double2 vr = sub2(x, y);
      const double rr = dot2(vr, vr);
      /* Evaluate kernel */
#if defined(CBIE)
      const double kern = (vr.x * dy.y - vr.y * dy.x) / rr;
#else
      //140225	const double kern = (nx.x * dy.y - nx.y * dy.x) / rr - 2.0 * (vr.x * nx.x + vr.y * nx.y) * (vr.x * dy.y - vr.y * dy.x) / (rr * rr);
      const double kern = ((nx.x * dy.y - nx.y * dy.x) - 2.0 * (vr.x * nx.x + vr.y * nx.y) * (vr.x * dy.y - vr.y * dy.x) / rr) / rr;
#endif
      /* Compute the potential at y */
      double uval = 0.0;
      for (int l = 0; l <= p; l ++) { // then, k-p<=k-l<=k, i.e., j<=k-l<=j+p
	const double bsp = closedcurve_basis(p, n, t, k - l, tnow);
	uval += bsp * xvec[joff + (k - l) % (n - p)];
      }
      sum += gw[ig] * Jg * kern * uval;
    } // ig
    
#endif
    
    //140325    *dval += sum;            
    *dval = sum;

  } // integral
  
  
  //140325  } // j
  
  *dval *= PI2I;
  
#ifdef MYINTDE
  myintde_free(intde);
#endif
}



void matvec_nbie(const GaussOne *g, const int ncc, const Closedcurve *cc, const int ncol, const double *tc, const double2 *xc, const double *freeterm, const double *xvec, double *yvec)
{
  /* Execute the matrix-vector product in terms of Neumann problem,
     where the normal derivative of potential is assumed to be zero on
     the boundary */

  /* Initialise */  
  for (int ix = 0; ix < ncol; ix ++) {
    yvec[ix] = 0.0;
  }
  
  int ioff = 0; // offset for collocation points
  
  /* Loop over closed curves */
  for (int icc = 0; icc < ncc; icc ++) {

    const int pi = cc[icc].p; // degree
    const int ni = cc[icc].n; // nuxmber of control points
    const double *ti = cc[icc].t; // knot values
    const double2 *cpi = cc[icc].cp; // control points

    /* Loop over collocation points */
#if defined(_OPENMP)
#pragma omp parallel for // OpenMP DEFINED LOOP WAS PARALLELIZED.
#endif
    for (int i = 0; i < ni - pi; i ++) {
      
      const int ix = ioff + i; // 0<=ix<ncol
      const double2 x = xc[ix]; // collocation point
      const double2 dx = closedcurve_position_derivative(pi, ni, ti, cpi, tc[ix]);
      const double2 nx = scale2(1.0 / norm2(dx), set2(dx.y, - dx.x)); // outward normal vector at x
      
      int joff = 0; // offset for igogeometric elements
      
      /* Loop over closed curves */
      for (int jcc = 0; jcc < ncc; jcc ++) {
	
	const int p = cc[jcc].p; // degree
	const int n = cc[jcc].n; // number of control points
	const double *t = cc[jcc].t; // knot values
	const double2 *cp = cc[jcc].cp; // control points
		
	//140325	double dval;
	//140325	isolap2d_direct_nbie(g, p, n, t, cp, tc, xc, icc, jcc, &dval, xvec, ioff, joff, i, pi, ni, ti, cpi);
	//140325		
	//140325	/* Accumulation */
	//140325	yvec[ix] = dval;

	/* Loop over isogeometric elements, i.e. section I_j=[t_{j+p},t_{j+p+1}) */
	for (int j = 0; j < n - p; j ++) {

	  double dval;
	  isolap2d_direct_nbie(g, p, n, t, cp, icc, jcc, &dval, xvec, ioff, joff, i, x, nx, j);
	  yvec[ix] += dval;

	} // j
	
	joff += n - p;
	
      } // jcc

#if defined(CBIE)
      /* Add free term */
      yvec[ix] += freeterm[ix] * closedcurve_potential(pi, ni, ti, &(xvec[ioff]), tc[ix]);
#endif

    } // i

    ioff += ni - pi;

  } // icc

}


//140318void isolap2d_matvec_nbie(const int N, const void *A, const double *xvec, double *yvec)
//140318{ 
//140318  Isolap2d *p = (Isolap2d *)A;
//140318  
//140318  matvec_nbie(p->gauss, p->ncc, p->cc, N, p->tc, p->xc, p->freeterm, xvec, yvec);
//140318}



//140325void isolap2d_direct_nbie_fmm(const GaussOne *g, const int p, const int n, const double *t, const double2 *cp, const double *tc, const double2 *xc, const int icc, const int jcc, double *dval, const double *xvec, const int ioff, const int joff, const int i, const int j, const int pi, const int ni, const double *ti, const double2 *cpi)
//140325{  
//140325  /* Global index of collocation point */
//140325  const int ix = ioff + i; // i is local
//140325  
//140325  /* Collocation point x and its derivative dx/dt */
//140325  const double2 x = xc[ix];
//140325  const double2 dx = closedcurve_position_derivative(pi, ni, ti, cpi, tc[ix]);
//140325
//140325  /* Compute the outward normal vector at the collocation point */
//140325  double2 nx;
//140325  nx.x = dx.y / norm2(dx); // nx1
//140325  nx.y = - dx.x / norm2(dx); // nx2
//140325  
//140325  /* Set up for Guassian quadrature */
//140325  const int ng = g->n;
//140325  const double *gx = g->x;
//140325  const double *gw = g->w;
//140325  
//140325  /* Initialise */
//140325  *dval = 0.0;
//140325
//140325  //  /* Loop over sections; I_j=[t_{j+p},t_{j+p+1}) or [t_k,t_{k+1}), where 0<=j<n-p */
//140325  //  for (int j = 0; j < n - p; j ++) {
//140325
//140325  const int k = j + p; // p<=k<n
//140325
//140325   
//140325  /* Integral over I_j for x_i:=y(t_{i+p}) */
//140325  if (i == j && icc == jcc) { // Singular at the lefthand-side end of I_j or I_i
//140325    
//140325    const double x1 = 0.0; // lower bound of X
//140325    const double x2 = 1.0; // upper bound of X
//140325    
//140325    if (p == 2) {
//140325      
//140325      /* Compute vectors a and b (where k-1>=1 and k-2>=0) */
//140325      const double2 va = sub2(cp[k - 1], cp[k - 2]);
//140325      const double2 vb = sub2(add2(cp[k], cp[k - 2]), scale2(2.0, cp[k - 1]));
//140325      
//140325      const double u0 = xvec[joff + (k - 0) % (n - p)];
//140325      const double u1 = xvec[joff + (k - 1) % (n - p)];
//140325      const double u2 = xvec[joff + (k - 2) % (n - p)];
//140325      
//140325      const double U0 = (u1 + u2) / 2.0;
//140325      const double U1 = u1 - u2;
//140325      const double U2 = (u0 - 2.0 * u1 + u2) / 2.0;
//140325      
//140325      *dval += sekibun_p2(g, va, vb, U0, U1, U2, x1, x2);
//140325      
//140325    } else if (p == 3) {
//140325      
//140325      /* Compute vector a, b, and c (p=3); See Eq.(B.5) */
//140325      const double2 va = scale2(0.5, sub2(cp[k - 1], cp[k - 3]));
//140325      const double2 vb = sub2(add2(cp[k - 1], cp[k - 3]), scale2(2.0, cp[k - 2]));
//140325      const double2 vc = sub2(scale2(0.5, sub2(cp[k], cp[k - 3])), scale2(1.5, sub2(cp[k - 1], cp[k - 2])));       
//140325      
//140325      /* Compute potential (p=3); See Eq.(B.11) */      
//140325      const double u0 = xvec[joff + (k - 0) % (n - p)];
//140325      const double u1 = xvec[joff + (k - 1) % (n - p)];
//140325      const double u2 = xvec[joff + (k - 2) % (n - p)];
//140325      const double u3 = xvec[joff + (k - 3) % (n - p)];
//140325      
//140325      const double U0 = (u1 + 4.0 * u2 + u3) / 6.0;
//140325      const double U1 = (u1 - u3) / 2.0;
//140325      const double U2 = (u1 - 2.0 * u2 + u3) / 2.0;
//140325      const double U3 = (u0 - 3.0 * u1 + 3.0 * u2 - u3) / 6.0;
//140325      
//140325      *dval += sekibun_p3(g, va, vb, vc, U0, U1, U2, U3, x1, x2);
//140325      
//140325    } else if (p == 4) {
//140325      
//140325      const double2 va = scale2(1.0 / 24.0, add2(scale2(4.0, cp[k - 1]), add2(scale2(12.0, cp[k - 2]), add2(scale2(- 12.0, cp[k - 3]), scale2(- 4.0, cp[k - 4])))));
//140325      const double2 vb = scale2(1.0 / 24.0, add2(scale2(12.0, cp[k - 1]), add2(scale2(- 12.0, cp[k - 2]), add2(scale2(- 12.0, cp[k - 3]), scale2(12.0, cp[k - 4])))));
//140325      const double2 vc = scale2(1.0 / 24.0, add2(scale2(12.0, cp[k - 1]), add2(scale2(- 36.0, cp[k - 2]), add2(scale2(36.0, cp[k - 3]), scale2(- 12.0, cp[k - 4])))));
//140325      const double2 vd = scale2(1.0 / 24.0, add2(scale2(4.0, cp[k - 0]), add2(scale2(- 16.0, cp[k - 1]), add2(scale2(24.0, cp[k - 2]), add2(scale2(- 16.0, cp[k - 3]), scale2(4.0, cp[k - 4]))))));
//140325      
//140325      const double u0 = xvec[joff + (k - 0) % (n - p)];
//140325      const double u1 = xvec[joff + (k - 1) % (n - p)];
//140325      const double u2 = xvec[joff + (k - 2) % (n - p)];
//140325      const double u3 = xvec[joff + (k - 3) % (n - p)];
//140325      const double u4 = xvec[joff + (k - 4) % (n - p)];
//140325      
//140325      const double U0 = (           u1 + 11.0 * u2 + 11.0 * u3 +       u4) / 24.0;
//140325      const double U1 = (     4.0 * u1 + 12.0 * u2 - 12.0 * u3 - 4.0 * u4) / 24.0;
//140325      const double U2 = (     6.0 * u1 -  6.0 * u2 -  6.0 * u3 + 6.0 * u4) / 24.0;
//140325      const double U3 = (     4.0 * u1 - 12.0 * u2 + 12.0 * u3 - 4.0 * u4) / 24.0;
//140325      const double U4 = (u0 - 4.0 * u1 +  6.0 * u2 -  4.0 * u3 +       u4) / 24.0;
//140325      
//140325      *dval += sekibun_p4(g, va, vb, vc, vd, U0, U1, U2, U3, U4, x1, x2);
//140325      
//140325    } else { 
//140325      
//140325      /* p>=5 is not available */
//140325      abort();
//140325      
//140325    }
//140325    
//140325  } else if (i == (j + 1) % (n - p) && icc == jcc) { // Singular at the righthand-side of I_j or I_{i-1}=[t_{i+p-1},t_{i+p})
//140325    
//140325    const double x1 = - 1.0; // lower bound of X
//140325    const double x2 = 0.0; // upper bound of X
//140325    
//140325    if (p == 2) {
//140325      
//140325      /* Compute vectors a and b (where k-1>=1 and k-2>=0) */
//140325      const double2 va = sub2(cp[k], cp[k - 1]);
//140325      const double2 vb = sub2(add2(cp[k], cp[k - 2]), scale2(2.0, cp[k - 1]));
//140325      
//140325      const double u0 = xvec[joff + (k - 0) % (n - p)];
//140325      const double u1 = xvec[joff + (k - 1) % (n - p)];
//140325      const double u2 = xvec[joff + (k - 2) % (n - p)];
//140325      
//140325      const double U0 = (u0 + u1) / 2.0;
//140325      const double U1 = u0 - u1;
//140325      const double U2 = (u0 - 2.0 * u1 + u2) / 2.0;
//140325      
//140325      *dval += sekibun_p2(g, va, vb, U0, U1, U2, x1, x2);
//140325      
//140325    } else if (p == 3) { 
//140325      
//140325      /* Compute vector a, b, and c (p=3); See Eq.(B.23) */
//140325      const double2 va = scale2(0.5, sub2(cp[k], cp[k - 2]));
//140325      const double2 vb = sub2(add2(cp[k], cp[k - 2]), scale2(2.0, cp[k - 1]));
//140325      const double2 vc = sub2(scale2(0.5, sub2(cp[k], cp[k - 3])), scale2(1.5, sub2(cp[k - 1], cp[k - 2])));
//140325      
//140325      const double u0 = xvec[joff + (k - 0) % (n - p)];
//140325      const double u1 = xvec[joff + (k - 1) % (n - p)];
//140325      const double u2 = xvec[joff + (k - 2) % (n - p)];
//140325      const double u3 = xvec[joff + (k - 3) % (n - p)];
//140325      
//140325      const double U0 = (u0 + 4.0 * u1 + u2) / 6.0;
//140325      const double U1 = (u0 - u2) / 2.0;
//140325      const double U2 = (u0 - 2.0 * u1 + u2) / 2.0;
//140325      const double U3 = (u0 - 3.0 * u1 + 3.0 * u2 - u3) / 6.0;
//140325      
//140325      *dval += sekibun_p3(g, va, vb, vc, U0, U1, U2, U3, x1, x2);
//140325      
//140325    } else if (p == 4) {
//140325      
//140325      const double2 va = scale2(1.0 / 24.0, add2(scale2(4.0, cp[k - 0]), add2(scale2(12.0, cp[k - 1]), add2(scale2(- 12.0, cp[k - 2]), scale2(- 4.0, cp[k - 3])))));
//140325      const double2 vb = scale2(1.0 / 24.0, add2(scale2(12.0, cp[k - 0]), add2(scale2(- 12.0, cp[k - 1]), add2(scale2(- 12.0, cp[k - 2]), scale2(12.0, cp[k - 3])))));
//140325      const double2 vc = scale2(1.0 / 24.0, add2(scale2(12.0, cp[k - 0]), add2(scale2(- 36.0, cp[k - 1]), add2(scale2(36.0, cp[k - 2]), scale2(- 12.0, cp[k - 3])))));
//140325      const double2 vd = scale2(1.0 / 24.0, add2(scale2(4.0, cp[k - 0]), add2(scale2(- 16.0, cp[k - 1]), add2(scale2(24.0, cp[k - 2]), add2(scale2(- 16.0, cp[k - 3]), scale2(4.0, cp[k - 4]))))));
//140325      
//140325      const double u0 = xvec[joff + (k - 0) % (n - p)];
//140325      const double u1 = xvec[joff + (k - 1) % (n - p)];
//140325      const double u2 = xvec[joff + (k - 2) % (n - p)];
//140325      const double u3 = xvec[joff + (k - 3) % (n - p)];
//140325      const double u4 = xvec[joff + (k - 4) % (n - p)];
//140325
//140325      const double U0 = (      u0 + 11.0 * u1 + 11.0 * u2 +       u3) / 24.0;
//140325      const double U1 = (4.0 * u0 + 12.0 * u1 - 12.0 * u2 - 4.0 * u3) / 24.0;
//140325      const double U2 = (6.0 * u0 -  6.0 * u1 -  6.0 * u2 + 6.0 * u3) / 24.0;
//140325      const double U3 = (4.0 * u0 - 12.0 * u1 + 12.0 * u2 - 4.0 * u3) / 24.0;
//140325      const double U4 = (      u0 -  4.0 * u1 +  6.0 * u2 - 4.0 * u3 + u4) / 24.0;
//140325
//140325      *dval += sekibun_p4(g, va, vb, vc, vd, U0, U1, U2, U3, U4, x1, x2);
//140325
//140325    } else { 
//140325
//140325      /* p>=5 is not available */
//140325      abort();
//140325	
//140325    }
//140325      
//140325
//140325  } else { // Non-singular integral; See Eq.(2.23) in Sakabe's thesis (which looks to miss the factor '2' in the enumerator of the second term in the RHS)
//140325      
//140325
//140325    double sum = 0.0; // integral over I_j for x_i
//140325
//140325    /* Compute the ends of I_j */
//140325    const double tminus = t[k];
//140325    const double tplus = t[k + 1];
//140325
//140325    /* Jacobian wrt Gaussian quadrature (this is not Jacobian wrt the curve parameter) */
//140325    const double Jg = (tplus - tminus) / 2.0; // (t_{k+1}-t_{k})/2
//140325
//140325    for (int ig = 0; ig < ng; ig ++) {
//140325      /* Parameter corresponding to the ig th abscissa */
//140325      const double tnow = ((1.0 - gx[ig]) * tminus + (1.0 + gx[ig]) * tplus) / 2.0;
//140325      /* Integral point y and its derivative dy/dt */	      
//140325      const double2 y = closedcurve_position(p, n, t, cp, tnow);
//140325      const double2 dy = closedcurve_position_derivative(p, n, t, cp, tnow);
//140325      /* Compute r:=|x-y| */
//140325      const double2 vr = sub2(x, y);
//140325      const double rr = dot2(vr, vr);
//140325      /* Evaluate kernel */
//140325#if defined(CBIE)
//140325      const double kern = (vr.x * dy.y - vr.y * dy.x) / rr;
//140325#else
//140325      const double kern = ((nx.x * dy.y - nx.y * dy.x) - 2.0 * (vr.x * nx.x + vr.y * nx.y) * (vr.x * dy.y - vr.y * dy.x) / rr) / rr;
//140325#endif
//140325      /* Compute the sum for b-spline and potential */
//140325      double uval = 0.0;
//140325      for (int l = 0; l <= p; l ++) { // then, k-p<=k-l<=k, i.e., j<=k-l<=j+p
//140325	const double bsp = closedcurve_basis(p, n, t, k - l, tnow);
//140325	uval += bsp * xvec[joff + (k - l) % (n - p)];
//140325      }
//140325      sum += gw[ig] * Jg * kern * uval;
//140325    } // ig
//140325
//140325    *dval += sum;            
//140325
//140325  } // integral
//140325    
//140325  //140318  } // j
//140325  
//140325  *dval *= PI2I;
//140325  
//140325}
