#include "options.h"


void isolap2d_options(int argc, char **argv, char **dxffile, char **infile)
{
  /* Check all the options and arguments */

  {
    int i = 0;
    INFO("%s", argv[i ++]);
    while (i < argc) {
      fprintf(stderr, " %s", argv[i ++]);
    }
    fprintf(stderr, "\n");
  }

  /* Set defaults */

  *dxffile = NULL;
  *infile = NULL;

  /* Check options */

  int c;
  FILE *fp;

  while ((c = getopt(argc, argv, "d:i:")) != - 1) {
    switch (c) {

    case 'd':

      if (optarg == 0) {
	MESG("Speficify a DXF file. Exit.\n");
	exit(EXIT_FAILURE);
      } else {	
	*dxffile = optarg;
	INFO("dxffile=%s\n", *dxffile);
      }

      break;

    case 'i':

      if (optarg == 0) {
	MESG("Speficify an infile. Exit.\n");
	exit(EXIT_FAILURE);
      } else {	
	*infile = optarg;
	INFO("infile=%s\n", *infile);

	fp = fopen(*infile, "r");
	if (fp == NULL) {
	  INFO("Fail to open %s. Exit.\n", *infile);
	  exit(EXIT_FAILURE);
	} else {
	  INFO("Pass to opening test of %s.\n", *infile);
	}
	  
      }

      break;

    default:

      MESG("Invalid option. Exit.\n");
      exit(EXIT_FAILURE);

    }
  }

}
