#include "elapsed.h"

#ifndef ENABLE_GETTIMEOFDAY

#include <stdio.h>
#include <stdlib.h>

#ifdef __i386__
#define RDTSC(X) asm volatile ("rdtsc" : "=A" (X))
#elif __x86_64__
#define RDTSC(X) asm volatile ("rdtsc; shlq $32, %%rdx; orq %%rdx, %%rax" : "=a" (X) :: "%rdx")
#else
#error Unknown architecture. RDTSC cannot be defined.
#endif

#if !defined(CPU_CLOCK_GHZ)
#error CPU_CLOCK_GHZ is undefined.
#endif

double elapsed(void)
{
  unsigned long long t;
  RDTSC(t);
  return (double)t / CPU_CLOCK_GHZ * 1.0e-9;
}

#else

/* http://wwweic.eri.u-tokyo.ac.jp/computer/manual/lx/prog/prog3.html */
#include <sys/time.h>
#include <stdio.h>

double elapsed(void)
{
  struct timeval tp;
  struct timezone tzp;
  gettimeofday(&tp,&tzp);
  return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

#endif

double elapsed_(void)
{
  return elapsed();
}

