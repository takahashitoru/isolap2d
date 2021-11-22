#ifndef ENVS_H
#define ENVS_H

#include <stdio.h>
#include <stdlib.h>

#ifndef GETENV
#define GETENV(msg) if (getenv(msg) != NULL) fprintf(stderr, "# %s: %s = %s\n", __FUNCTION__, msg, getenv(msg))
#endif

#ifdef __cplusplus
extern "C" {
#endif
  void envs(void);
#ifdef __cplusplus
}
#endif

#endif /* ENVS_H */
