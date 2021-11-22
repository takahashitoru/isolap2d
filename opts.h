#ifndef OPTS_H
#define OPTS_H

#include <stdio.h>
#include <stdlib.h>

#include "debugmacros.h"

#ifndef MESSAGE
#define MESSAGE(msg) fprintf(stderr, "# %s: %s\n", __FUNCTION__, msg)
#endif
#ifndef ERRMESG
#define ERRMESG(msg) fprintf(stderr, "# %s: %s Exit.\n", __FUNCTION__, msg); exit(1)
#endif
#ifndef PRINT
#define PRINT(fmt, ...) fprintf(stderr, "# %s: " fmt "\n", __FUNCTION__, __VA_ARGS__)
#endif
#ifndef SOP
#define SOP(option, value) INFO("%s = %s\n", option, value)
#endif
#ifndef IOP
#define IOP(option, value) INFO("%s = %d\n", option, value)
#endif
#ifndef DOP
#define DOP(option, value) INFO("%s = %le\n", option, value)
#endif

#ifdef __cplusplus
extern "C" {
#endif
  void opts(void);
#ifdef __cplusplus
}
#endif

#endif /* OPTS_H */
