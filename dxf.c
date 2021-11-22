#include "dxf.h"

/* Inputed .dxf is assumed to be created as "closed splines" of degree
   2 or more by LibrCAD and contains one or more B-spline closed
   curves.

   (Input) dxffile: File in the DXF format.

   (Output) p: Degree of B-splines; two or more.

   (Output) n: Number of control points, which is equal to
   "num_control_points + degree" because of wrapping control points.
   
   (Output) cp: Coordinates of n control points.
*/

void dxf_read_bspline(const char *dxffile, const int maxbsp, int *nbsp, Bspline *bsp)
{

  FILE *fptr = fopen(dxffile, "r");
  if (fptr == NULL) {
    INFO("Fail to open %s. Exit.\n", dxffile);
    exit(EXIT_FAILURE);    
  }

  int ibsp = - 1; // initialise the number of B-splines (not 0, but -1)
  char *code = (char *)malloc(DXF_MAX_COLUMNS * sizeof(char));
  char *value = (char *)malloc(DXF_MAX_COLUMNS * sizeof(char));
  int c;
  int entity_spline = DXF_FALSE;
  
  do {

    fgets(code, DXF_MAX_COLUMNS, fptr) ; // read a code
    c = atoi(code);
    fgets(value, DXF_MAX_COLUMNS, fptr); // read its value
    
    if (c == 0) {

      if (strncmp(value, "SPLINE", 6) == 0) { // begin SPLINE entity

	ibsp ++; // set the counter of B-splines firstly
	entity_spline = DXF_TRUE;

	if (ibsp >= maxbsp) {
	  INFO("maxpsp = %d is too small.\n", maxbsp);
	  exit(EXIT_FAILURE);
	}

	/////////////////////////////////////////////////////////////////////////
	DBG("%d: ibsp=%d: Enter SPLINE\n", c, ibsp);
	/////////////////////////////////////////////////////////////////////////
	
      } else {

	entity_spline = DXF_FALSE;

      }
      
    } else if (entity_spline == DXF_TRUE && c == 70) { // Spline flag
      
      int flag = atoi(value);
      
      if (flag != 11) {
	INFO("%d: ibsp=%d: Spline flag looks strange: %d\n", c, ibsp, flag);
	exit(EXIT_FAILURE);
      } 

      bsp[ibsp].flag = flag;

      /////////////////////////////////////////////////////////////////////////
      DBG("%d: ibsp=%d: flag=%d\n", c, ibsp, flag);
      /////////////////////////////////////////////////////////////////////////
      
    } else if (entity_spline == DXF_TRUE && c == 71) { // Degree of the spline curve
      
      int degree = atoi(value);
      
      if (degree < 2) {
	INFO("%d: ibsp=%d: Degree of the spline curve looks invalid: %d\n", c, ibsp, degree);
	exit(EXIT_FAILURE);
      } 
      
      bsp[ibsp].degree = degree;
      
      /////////////////////////////////////////////////////////////////////////
      DBG("%d: ibsp=%d: degree=%d\n", c, ibsp, degree);
      /////////////////////////////////////////////////////////////////////////

    } else if (entity_spline == DXF_TRUE && c == 72) { // Number of knots
      
      /////////////////////////////////////////////////////////////////////////
      //      DBG("%d: ibsp=%d: Number of knots is ignored\n", c, ibsp);
      /////////////////////////////////////////////////////////////////////////
      
    } else if (entity_spline == DXF_TRUE && c == 73) { // Number of control points
      
      int num_control_points = atoi(value);
      
      if (num_control_points <= 0) {
	INFO("%d: ibsp=%d: Number of control points looks invalid: %d\n", c, ibsp, num_control_points);
	exit(EXIT_FAILURE);
      }

      /////////////////////////////////////////////////////////////////////////
      DBG("%d: ibsp=%d: num_cotrol_points=%d\n", c, ibsp, num_control_points);
      /////////////////////////////////////////////////////////////////////////

      bsp[ibsp].m = num_control_points;
      bsp[ibsp].x = (double *)malloc(num_control_points * sizeof(double));
      bsp[ibsp].y = (double *)malloc(num_control_points * sizeof(double));
      bsp[ibsp].z = (double *)malloc(num_control_points * sizeof(double)); // always zero for 2D
      bsp[ibsp].cx = 0; // initialise the counter
      bsp[ibsp].cy = 0; // initialise the counter
      bsp[ibsp].cz = 0; // initialise the counter

    } else if (entity_spline == DXF_TRUE && c == 40) { // Knot value

      /////////////////////////////////////////////////////////////
      //      DBG("%d: ibsp=%d: Knot value is ignored.\n", c, ibsp);
      /////////////////////////////////////////////////////////////

    } else if (entity_spline == DXF_TRUE && c == 10) { // X value of control point
      (bsp[ibsp].x)[(bsp[ibsp].cx) ++] = atof(value);
      
    } else if (entity_spline == DXF_TRUE && c == 20) { // Y value of control point
      (bsp[ibsp].y)[(bsp[ibsp].cy) ++] = atof(value);
      
    } else if (entity_spline == DXF_TRUE && c == 30) { // Z value of control point
      (bsp[ibsp].z)[(bsp[ibsp].cz) ++] = atof(value);
      
    } else {
      ////////////////////////////////////////////////////////
      //      DBG("code:  %s", code);
      //      DBG("value: %s", value);
      ////////////////////////////////////////////////////////
    }
      
  } while (c != 0 || strncmp(value, "EOF", 3) != 0); // exit if the end of file
  
  fclose(fptr);

  *nbsp = ibsp + 1; // number of B-splines read 

  if (*nbsp == 0) {
    INFO("Fail to read any SPLINE from %s. Exit.\n", dxffile);
    exit(EXIT_FAILURE);
  } else {
    INFO("Read %d SPLINE(s) from %s.\n", *nbsp, dxffile);
  }
  
  ///////////////////////////////////////////////////////////////////////
#if(0)
  for (int ibsp = 0; ibsp < *nbsp ; ibsp ++) {
    int p = bsp[ibsp].degree;
    int m = bsp[ibsp].m;
    DBG("ibsp=%d: p=%d m=%d\n", ibsp, p, m);
    for (int i = 0; i < m; i ++) {
      DBG("i=%d: x=%f y=%f z=%f\n", i, bsp[ibsp].x[i], bsp[ibsp].y[i], bsp[ibsp].z[i]);
    }
  }
#endif
  ///////////////////////////////////////////////////////////////////////
    
  free(code);
  free(value);
}

void dxf_free_bspline(const int nbsp, Bspline *bsp)
{
  for (int ibsp = 0; ibsp < nbsp; ibsp ++) {
    free(bsp[ibsp].x);
    free(bsp[ibsp].y);
    free(bsp[ibsp].z);
  }
  free(bsp);
}
