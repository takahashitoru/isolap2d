#include "bincoef.h"
//double bincoef(int n, int k)
//{
//  /*
//    Compute binomial coefficient: 
//    http://en.wikipedia.org/wiki/Binomial_coefficient
//  */
//
//  assert(n >=0 && 0 <= k && k <= n);
//  
//  if (k > n - k) { // take advantage of symmetry
//    k = n - k;
//  }
//  double c = 1;
//  for (int i = 0; i < k; i ++) {
//    c = c * (n - (k - (i + 1)));
//    c = c / (i + 1);
//  }
//  return c;
//}


double *bincoef_malloc(const int m)
{
  return (double *)malloc(m * m * sizeof(double));
}

void bincoef_free(double *b)
{
  free(b);
}

void bincoef(const int m, double *b)
{
  /*
    Compute binomial coefficients b[m * n + k]=(n,k), where 0<=n<m and 0<=k<m.
    Allocate like double *b = (double *)malloc(m * m * sizeof(double));
    See http://en.wikipedia.org/wiki/Binomial_coefficient
  */

  for (int n = 0; n < m; n ++) {
    b[m * n + 0] = 1.0; // (n,0)=1 for n>=0
  }    
  for (int k = 1; k < m; k ++) {
    b[m * 0 + k] = 0.0; // (0,k)=0 for k>0
  }    
  for (int n = 1; n < m; n ++) {
    for (int k = 1; k < m; k ++) { // not k<=m, but k<m
      b[m * n + k] = b[m * (n - 1) + (k - 1)] + b[m * (n - 1) + k];
    }
  }
}
