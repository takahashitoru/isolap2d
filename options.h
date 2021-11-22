#ifndef OPTIONS_H
#define OPTIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include "debugmacros.h"

#ifdef __cplusplus
extern "C" {
#endif

  void isolap2d_options(int argc, char **argv, char **dxffile, char **infile);

#ifdef __cplusplus
}
#endif

#endif /* OPTIONS_H */
