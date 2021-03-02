/***************************************************************
 * This file contains procedure performingNewton's polynomial
 * interpolation. 
 */
#include <stdio.h>
#include <stdlib.h>


double  *c=NULL;


/*---------------------------------------------------------
 * Procedure to compute the coefficients of Newton
 * polynomial.
 *   data[][2]: sample points, (xi, yi)
 *   n: largest index of sample point.
 * Return the interpolated value.
 */
void comp_Newton_coef(double data[][2], int n)
{
	int     i, j;


	// Create space for the coef's.
	c = (double *) malloc(sizeof(double)*(n+1));

	//initialize the coef's.
	for(i=0;i<=n;i++) c[i] = data[i][1];

	//Recursively compute the coefficients.
	for(i=0;i<=n;i++){
		for(j=n;j>i;j--)  // Modify c_j = (c_j-c_i)/(xj-xi)
			c[j] = (c[j]-c[i])/(data[j][0]-data[i][0]);
	}
    for(i=0;i<=n;i++)
		fprintf(stderr," c%d = %lf\n", i, c[i]);
}


/*-------------------------------------------------------------
 * Procedure to perform Newton polynomial interpolation.
 *   data[][2]: the sample points,
 *   t: the input parameter value,
 *   n: last index of sample pnt.
 * Return the interpolayed value.
 */
double Newton_interpolate(double data[][2], double t, int n)
{
	double   sum=0.0;
    int      i;


	//Compute polynomial value by using Horner's algm.
	sum = c[n];
	for(i=n-1;i>=0;i--)
		sum = sum*(t-data[i][0]) + c[i];
	return(sum);
}
