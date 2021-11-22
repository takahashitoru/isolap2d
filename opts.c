#include "opts.h"

void opts(void)
{
#ifdef __DATE__
  SOP("__DATE__", __DATE__);
#endif
#ifdef __TIME__
  SOP("__TIME__", __TIME__);
#endif

#ifdef __linux__
  MESSAGE("__linux__");
#endif

#ifdef __i386__
  MESSAGE("__i386__");
#endif
#ifdef __x86_64__
  MESSAGE("__x86_64__");
#endif

#ifdef __SSE__
  MESSAGE("__SSE__");
#endif
#ifdef __SSE2__
  MESSAGE("__SSE2__");
#endif
#ifdef __SSE3__
  MESSAGE("__SSE3__");
#endif
#ifdef __SSSE3__
  MESSAGE("__SSSE3__");
#endif
#ifdef __SSE4_1__
  MESSAGE("__SSE4_1__");
#endif
#ifdef __SSE4_2__
  MESSAGE("__SSE4_2__");
#endif
#ifdef __AVX__
  MESSAGE("__AVX__");
#endif

#ifdef __ICC
  IOP("__ICC", __ICC);
#endif
#ifdef __INTEL_COMPILER_BUILD_DATE
  IOP("__INTEL_COMPILER_BUILD_DATE", __INTEL_COMPILER_BUILD_DATE);
#endif

#ifdef __GNUC__
  IOP("__GNUC__", __GNUC__);
#endif
#ifdef __GNUC_MINOR__
  IOP("__GNUC_MINOR__", __GNUC_MINOR__);
#endif
#ifdef __GNUC_PATCHLEVEL__
  IOP("__GNUC_PATCHLEVEL__", __GNUC_PATCHLEVEL__);
#endif

#ifdef _OPENMP
  MESSAGE("_OPENMP");
#endif

#ifdef MYDEBUG
  MESSAGE("MYDEBUG");
#endif

#ifdef CIRCLE
  MESSAGE("CIRCLE");
#endif

#ifdef FMM
  MESSAGE("FMM");
#endif

#ifdef FUKUI
  MESSAGE("FUKUI");
#endif

#ifdef HEXA
  MESSAGE("HEXA");
#endif

#ifdef SET_FREE_TERM
  MESSAGE("SET_FREE_TERM");
#endif

#ifdef MYINTDE
  MESSAGE("MYINTDE");
#ifdef MYINTDE_EPS
  DOP("MYINTDE_EPS", MYINTDE_EPS);
#endif
#ifdef MYINTDE_LENAW
  IOP("MYINTDE_LENAW", MYINTDE_LENAW);
#endif
#endif

#ifdef ADAPTIVE_GAUSS
  MESSAGE("ADAPTIVE_GAUSS");
#ifdef ADAPTIVE_GAUSS_NMIN
  IOP("ADAPTIVE_GAUSS_NMIN", ADAPTIVE_GAUSS_NMIN);
#endif
#ifdef ADAPTIVE_GAUSS_NMAX
  IOP("ADAPTIVE_GAUSS_NMAX", ADAPTIVE_GAUSS_NMAX);
#endif
#ifdef ADAPTIVE_GAUSS_EPS
  DOP("ADAPTIVE_GAUSS_EPS", ADAPTIVE_GAUSS_EPS);
#endif
#endif

}
