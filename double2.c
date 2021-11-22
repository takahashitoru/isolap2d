#include <double2.h>

double dot2(const double2 u, const double2 v)
{
  return u.x * v.x + u.y * v.y;
}

double norm2(const double2 u)
{
  return sqrt(dot2(u, u));
}

double2 add2(const double2 u, const double2 v)
{
  double2 w;
  w.x = u.x + v.x;
  w.y = u.y + v.y;
  return w;
}

double2 sub2(const double2 u, const double2 v)
{
  double2 w;
  w.x = u.x - v.x;
  w.y = u.y - v.y;
  return w;
}

double2 scale2(const double a, const double2 u)
{
  double2 w;
  w.x = a * u.x;
  w.y = a * u.y;
  return w;
}

double vdot2(const double2 u, const double2 v)
{
  return u.x * v.y - u.y * v.x;
}

double2 set2(const double x, const double y)
{
  double2 d;
  d.x = x;
  d.y = y;
  return d;
}
