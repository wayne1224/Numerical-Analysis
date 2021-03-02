/***************************************************************
 * This file contains procedure performing Lagrange's polynomial
 * interpolation. 
 */
#include <stdio.h>
#include <stdlib.h>


/*---------------------------------------------------------
 * Procedure to compute the integral by using Simpson's rule.
 *   data[][2]: sample points, (xi, yi)
 *   h: gap between sample points,
 *   n: largest index of sample point.
 * n must be an even number.
 * Return the integral value.
 */
double simpson(double data[][2], double h, int n)
{
	int     i;
	double  sum=0.0;


	sum = data[0][1] + data[n][1];  // sum = f0 + fn;

	// sum += 2*fi, 1<= i <= n-1
	for(i=1;i<n-1;i+=2)
		sum = sum + 4.0*data[i][1] + 2.0*data[i+1][1];
	sum += 4.0*data[n-1][1];

    // scale integral by h/3.
    sum = sum*h/3.0;
	return(sum);
}

