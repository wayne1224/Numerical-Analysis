/**************************************************************
 * This file contains definitions of onstants, procedure headers and other information
 *  related to the programs.
 */

#define   MY_PI     3.141592653

#define   RADIUS   10.0

//Number of interpolation points
#define   NUM_INTER_PNTS   1000



/##########################################################
 # Procedure headers
 # */


/*-----------------------------------------------------------------------------------------
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
void  gen_circle_pnts(double *t, double  *x, double *y, int N, double r);

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
void Lagrange(int N, double *t, double *f, double *p);

