/*
DE-Quadrature
Numerical Automatic Integrator for Improper Integral
    method    : Double Exponential (DE) Transformation
    dimension : one
    table     : use
functions
    intde  : integrator of f(x) over (a,b).
    intdei : integrator of f(x) over (a,infinity), 
                 f(x) is non oscillatory function.
    intdeo : integrator of f(x) over (a,infinity), 
                 f(x) is oscillatory function.
*/

/*
intde
    [description]
        I = integral of f(x) over (a,b)
    [declaration]
        void intdeini(int lenaw, double tiny, double eps, double *aw);
        void intde(double (*f)(double), double a, double b, 
            double *aw, double *i, double *err);
    [usage]
        intdeini(lenaw, tiny, eps, aw);  // initialization of aw
        ...
        intde(f, a, b, aw, &i, &err);
    [parameters]
        lenaw     : length of aw (int)
        tiny      : minimum value that 1/tiny does not 
                    overflow (double)
        eps       : relative error requested (double)
        aw        : points and weights of the quadrature 
                    formula, aw[0...lenaw-1] (double *)
        f         : integrand f(x) (double (*f)(double))
        a         : lower limit of integration (double)
        b         : upper limit of integration (double)
        i         : approximation to the integral (double *)
        err       : estimate of the absolute error (double *)
    [remarks]
        initial parameters
            lenaw > 1000, 
            IEEE double :
                lenaw = 8000;
                tiny = 1.0e-307;
        function
            f(x) needs to be analytic over (a,b).
        relative error
            eps is relative error requested excluding 
            cancellation of significant digits.
            i.e. eps means : (absolute error) / 
                             (integral_a^b |f(x)| dx).
            eps does not mean : (absolute error) / I.
        error message
            err >= 0 : normal termination.
            err < 0  : abnormal termination.
                       i.e. convergent error is detected :
                           1. f(x) or (d/dx)^n f(x) has 
                              discontinuous points or sharp 
                              peaks over (a,b).
                              you must divide the interval 
                              (a,b) at this points.
                           2. relative error of f(x) is 
                              greater than eps.
                           3. f(x) has oscillatory factor 
                              and frequency of the oscillation 
                              is very high.
*/


//#include <math.h>

#include "myintde2.h"

void intdeini(int lenaw, double tiny, double eps, double *aw)
{
    /* ---- adjustable parameter ---- */
    double efs = 0.1, hoff = 8.5;
    /* ------------------------------ */
    int noff, nk, k, j;
    double pi2, tinyln, epsln, h0, ehp, ehm, h, t, ep, em, xw, wg;
    
    pi2 = 2 * atan(1.0);
    tinyln = -log(tiny);
    epsln = 1 - log(efs * eps);
    h0 = hoff / epsln;
    ehp = exp(h0);
    ehm = 1 / ehp;
    aw[2] = eps;
    aw[3] = exp(-ehm * epsln);
    aw[4] = sqrt(efs * eps);
    noff = 5;
    aw[noff] = 0.5;
    aw[noff + 1] = h0;
    aw[noff + 2] = pi2 * h0 * 0.5;
    h = 2;
    nk = 0;
    k = noff + 3;
    do {
        t = h * 0.5;
        do {
            em = exp(h0 * t);
            ep = pi2 * em;
            em = pi2 / em;
            j = k;
            do {
                xw = 1 / (1 + exp(ep - em));
                wg = xw * (1 - xw) * h0;
                aw[j] = xw;
                aw[j + 1] = wg * 4;
                aw[j + 2] = wg * (ep + em);
                ep *= ehp;
                em *= ehm;
                j += 3;
            } while (ep < tinyln && j <= lenaw - 3);
            t += h;
            k += nk;
        } while (t < 1);
        h *= 0.5;
        if (nk == 0) {
            if (j > lenaw - 6) j -= 3;
            nk = j - noff;
            k += nk;
            aw[1] = nk;
        }
    } while (2 * k - noff - 3 <= lenaw);
    aw[0] = k - 3;
}


//void intde(double (*f)(double), double a, double b, double *aw, double *i, double *err)
static void intde(Intde *p, double (*f)(double, Intde *), double a, double b, double *aw, double *i, double *err)
{
    int noff, lenawm, nk, k, j, jtmp, jm, m, klim;
    double epsh, ba, ir, xa, fa, fb, errt, errh, errd, h, iback, irback;
    
    noff = 5;
    lenawm = (int) (aw[0] + 0.5);
    nk = (int) (aw[1] + 0.5);
    epsh = aw[4];
    ba = b - a;
    //    *i = (*f)((a + b) * aw[noff]);
    *i = (*f)((a + b) * aw[noff], p);
    ir = *i * aw[noff + 1];
    *i *= aw[noff + 2];
    *err = fabs(*i);
    k = nk + noff;
    j = noff;
    do {
        j += 3;
        xa = ba * aw[j];
	//        fa = (*f)(a + xa);
        fa = (*f)(a + xa, p);
	//        fb = (*f)(b - xa);
        fb = (*f)(b - xa, p);
        ir += (fa + fb) * aw[j + 1];
        fa *= aw[j + 2];
        fb *= aw[j + 2];
        *i += fa + fb;
        *err += fabs(fa) + fabs(fb);
    } while (aw[j] > epsh && j < k);
    errt = *err * aw[3];
    errh = *err * epsh;
    errd = 1 + 2 * errh;
    jtmp = j;
    while (fabs(fa) > errt && j < k) {
        j += 3;
	//        fa = (*f)(a + ba * aw[j]);
        fa = (*f)(a + ba * aw[j], p);
        ir += fa * aw[j + 1];
        fa *= aw[j + 2];
        *i += fa;
    }
    jm = j;
    j = jtmp;
    while (fabs(fb) > errt && j < k) {
        j += 3;
	//        fb = (*f)(b - ba * aw[j]);
        fb = (*f)(b - ba * aw[j], p);
        ir += fb * aw[j + 1];
        fb *= aw[j + 2];
        *i += fb;
    }
    if (j < jm) jm = j;
    jm -= noff + 3;
    h = 1;
    m = 1;
    klim = k + nk;
    while (errd > errh && klim <= lenawm) {
        iback = *i;
        irback = ir;
        do {
            jtmp = k + jm;
            for (j = k + 3; j <= jtmp; j += 3) {
                xa = ba * aw[j];
		//                fa = (*f)(a + xa);
                fa = (*f)(a + xa, p);
		//                fb = (*f)(b - xa);
                fb = (*f)(b - xa, p);
                ir += (fa + fb) * aw[j + 1];
                *i += (fa + fb) * aw[j + 2];
            }
            k += nk;
            j = jtmp;
            do {
                j += 3;
		//                fa = (*f)(a + ba * aw[j]);
                fa = (*f)(a + ba * aw[j], p);
                ir += fa * aw[j + 1];
                fa *= aw[j + 2];
                *i += fa;
            } while (fabs(fa) > errt && j < k);
            j = jtmp;
            do {
                j += 3;
		//                fb = (*f)(b - ba * aw[j]);
                fb = (*f)(b - ba * aw[j], p);
                ir += fb * aw[j + 1];
                fb *= aw[j + 2];
                *i += fb;
            } while (fabs(fb) > errt && j < k);
        } while (k < klim);
        errd = h * (fabs(*i - 2 * iback) + fabs(ir - 2 * irback));
        h *= 0.5;
        m *= 2;
        klim = 2 * klim - noff;
    }
    *i *= h * ba;
    if (errd > errh) {
        *err = -errd * (m * fabs(ba));
    } else {
        *err = *err * aw[2] * (m * fabs(ba));
    }
}



Intde *myintde_malloc()
{
  Intde *p = (Intde *)malloc(sizeof(Intde));
  p->lenaw = MYINTDE_LENAW;
  p->tiny = MYINTDE_TINY;
  p->eps = MYINTDE_EPS;
  return p;
}

void myintde_init(Intde *p)
{
  p->aw = (double *)malloc(p->lenaw * sizeof(double));
  intdeini(p->lenaw, p->tiny, p->eps, p->aw);
}


/*
  f         : integrand f(x) (double (*f)(double))
  a         : lower limit of integration (double)
  b         : upper limit of integration (double)
  i         : approximation to the integral (double *)
  err       : estimate of the absolute error (double *)
*/

void myintde(Intde *p, double (*f)(double, Intde *), double a, double b, double *i, double *err)
{
  intde(p, f, a, b, p->aw, i, err);
}


void myintde_free(Intde *p)
{
  free(p->aw);
  free(p);
}
