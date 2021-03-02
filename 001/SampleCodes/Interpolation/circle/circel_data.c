#include <stdio.h>
#include <math.h>

#include "mydefinition.h"

/************************************************************
 * Procedure to generate N sample points uniformly distributed in a circle.
 * The circel is entered at (0, 0)
 * Input:
 *    N: number of samples,
 *    r: radius of the circel.
 * Output:
 *    *t: array of parametric values (angles in radiens),
 *    *x, *y: arrays of x- & y-coordinates.
 *  Called by main() in main.c.
 */
void  gen_circle_pnts(double *t, double  *x, double *y, int N, double r)
{
   int  i;
   double   dTheta, theta;

   //compute the incremental value of angle in radians.
   dTheta = 2.0*MY_PI/((double)N);

  //Generate the x- and y-coordinates of the sample points.
   for(i=0;i<N;i++){
      theta = i*dTheta;
      t[i] = theta;
      x[i] = r*cos(theta);
      y[i] = r*sin(theta);
   }
    return;
} 