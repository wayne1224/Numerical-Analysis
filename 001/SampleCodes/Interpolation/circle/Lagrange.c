/*****************************************************************************
 * This file contains a Lagrange interpolation procedure which is used to generate interpolation
 *  points of the circel.
 */
#include <stdio.h>
#include "mydefinition.h"



/*-----------------------------------------------------------------------------------------------------
 * This procedure produces M interpolation points from input samples.
 *  Input:
 *      N: number of samples,
 *     *t:   array of parametric values,
 *     *f:   array of function values.
 * Output:
 *     *p: array of interpolations.
 * The number of interpolations = NUM_INTER_PNTS, defined in "mydefinition.h".
 */
void Lagrange(int N, double *t, double *f, double *p)
{
   int   i, j;
   double   term, sum, tmin, tmax, dt, tVal;

   //Compute the incremental value of the parameter t.
   tmax = t[N-1];
   tmin = t[0];
   dt = (tmax-tmin)/N;
   //Generate the interpolation points.
   for(k=0;k<NUM_INTER_PNTS;k++){ //for eah interpolation point
      tVal = k*dt;  //The parametric value of the interpolation point
      sum = 0.0;
      for(i=0;i<N;i++){ //for each term
         term = 1.0;
         for(j=0;j<N;j++) //Calculate the term
            if(i!=j) term = term*(tVal-t[j])/(t[i]-t[j]);
         sum = sum + term; //Acumulate the term
      }//end_for(i)
      p[k] = sum;
   }//end_for(k)
}

