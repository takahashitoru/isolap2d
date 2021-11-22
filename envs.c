#include "envs.h"

void envs(void)
{
  GETENV("PWD");
  GETENV("HOSTNAME");
#if defined(_OPENMP)
  GETENV("OMP_NUM_THREADS");
#endif
}

