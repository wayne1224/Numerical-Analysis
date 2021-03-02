/****************************************************************
 * This file contains procedure and data structures for computing
 * eiegn-values of symmetric tri-diagonal mtx which are obtained
 * from the Householder's transformation or Lanzcos method.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*--- Define the intervals for searching the eigen-values. ---*/
typedef struct{
	double y, z;
} myinterval_t;

myinterval_t  interval[1000]; // The intervals.
int         numInterval=0;  // The num. of intervals.
double      ymin, zmax;  // The maximum range of the interval.

// The values of the characteristic polynomial.
double  *p;



void comp_character_poly(double **T, int n, double x)
{
	int   i;

	p[0] = 1.0;
	p[1] = T[0][0] - x;
	for(i=2;i<=n;i++){
		p[i] = (T[i-1][i-1]-x)*p[i-1] - T[i-1][i-2]*T[i-1][i-2]*p[i-2];
	}
}

/*---------------------------------------------------------------
 * Procedure to evaluate the number of eigen-values less than a
 * threshold.
 * Input:
 *   T, n: the matrix and its dimension,
 *   x: the threshold.
 * Return:
 *   the number of eigen-values.
 */
int comp_num_eigen(double **T, int n, double x)
{
	int  i;
	int  change=0;

  // Compute the value of the characteristic polynomial.
	comp_character_poly(T, n, x);
	for(i=1;i<=n;i++){
		if(p[i-1]>=0 && p[i]<0) change ++;
		else if(p[i-1]<0 && p[i]>=0) change ++;
	}
	return (change);
}


/*---------------------------------------------------------------
 * Procedure to compute the supreme interval of eiegen values.
 * Input:
 *   T: a tri-diag mtx,
 *   n: dimension of T.
 * Output:
 *   y, z: y<z, the interval=[y,z].
 */
void comp_max_interval(double **T, int n, double *y, double *z)
{
	double   zmax, ymin, t0, t1;
	int      i;

	zmax = -1.0e7;
	ymin = 1.0e7;
	for(i=0;i<n;i++){
		if(i>0 && i<n-1){
			t0 = T[i][i] - fabs(T[i][i-1]) - fabs(T[i][i+1]);
			t1 = T[i][i] + fabs(T[i][i-1]) + fabs(T[i][i+1]);
		}else if(i==0){
			t0 = T[0][0] - fabs(T[0][1]);
			t1 = T[0][0] + fabs(T[0][1]);
		}else if(i==n-1){
			t0 = T[i][i] - fabs(T[i][i-1]);
			t1 = T[i][i] + fabs(T[i][i-1]);
		}
		if(t0<ymin) ymin = t0;
		if(t1>zmax) zmax = t1;
	}
	*y = ymin;
	*z = zmax;
	fprintf(stderr," The supreme interval: %lf - %lf\n", ymin, zmax);
}


/*--------------------------------------------------------------------
 * Procedure to output an eigen-value interval.
 *  y, z: the 2 ends of the interval.
 *  idx: index of the interval.
 */
void output_interval(double y, double z, int idx)
{
	interval[idx].y = y;
	interval[idx].z = z;
}


/*-----------------------------------------------------------------------
 * Recursive procedure to compute eigen-value intervals.
 * Input:
 *   numY, numZ: #(eigen-vals) < y and z.
 *   y, z: interval ends.
 */
void split_interval(double **T, int n, int numY, int numZ, double y, double z)
{
	double   x;
	int      numX;

	if(numY>=numZ) return; // empty interval.
	if(numZ == numY +1){  // Exact interval with 1 eigen-val.
		output_interval(y, z, numZ-1);
		return;
	}
	x = (y+z)/2.0;
	numX = comp_num_eigen(T, n, x);
	split_interval(T, n, numY, numX, y, x);
	split_interval(T, n, numX, numZ, x, z);
}


/*----------------------------------------------------------------------
 * Procedure to divide the supreme interval into small intervals.
 * Each small interval contains an eiegn-value.
 * Input:
 *   T: the tri-dia mtx,
 *   n: dim. of T.
 *   y, z: the 2 ends of the supreme interval.
 * Output:
 *   interval[]: the intervals
 *   numInterval: num. of intervals.
 */
void div_intervals(double **T, int n, double y, double z, myinterval_t *interval, int *numInterval)
{
	int  numY, numZ; // num. of eigenValues <=y and <= z.
	

    // Allocate space for the characteristic polynomial.
	p = (double *) malloc(sizeof(double)*(n+1));
	p[0] = 1.0;
	// Compute the number of eigen values less than y and z.
	numY = comp_num_eigen(T, n, y);
	numZ = comp_num_eigen(T, n, z);
	// Check the supreme interval and the resulted numbers.
	fprintf(stderr," There are %d and %d eigen-values less than %lf and %lf\n", numY, numZ, y, z);
	if(numY>=numZ){
		fprintf(stderr," The supreme interval is incorrect.\n");
		fprintf(stderr," Terminate the program.\n");
		return;
	}
	// Allocate spaces for the intervals.
//	interval = (myinterval_t *) malloc(sizeof(myinterval_t)*(numZ+1));
	*numInterval = numZ;

	// Recursively divide the intervals.
	split_interval(T, n, numY, numZ, y, z);
}



/*--------------------------------------------------------------------
 * Procedure to compute the eigenvalue by using the bisection method.
 * Input:
 *   T, n: the mtx and its dimension,
 *   y, z: the 2 ends of the interval of the eigenvalue.
 */
double bisect_method(double **T, int n, double y, double z)
{
	double   x, py, pz, px;
    
	
	comp_character_poly(T, n, y);
	py = p[n];
	comp_character_poly(T, n, z);
	pz = p[n];
	x = (y+z)/2.0;
    comp_character_poly(T, n, x);
	px = p[n];
	while((z-y)>0.000001){
		if(px==0.0) return x;
		if(px<0.0 && py<0.0){
			py = px;
			y = x;
		}else{
			pz = px;
			z = x;
		}
	    x = (y+z)/2.0;
        comp_character_poly(T, n, x);
	    px = p[n];
	}
	return(x);
}



/*---------------------------------------------------------------
 * The procedure to compute the eigen-values.
 * Input:
 *    T: a tri-diagonal matrix,
 *    n: dimension of the mtx.
 * Output:
 *    eigenVal[]: array of the eigen-values,
 *    numEigen: num. of eigen-values.
 */
void tri_diag_mtx_eigen(double **T, int n, double *eigenVal, int *numEigen)
{
	int  i;

	// Compute the max. interval of eigen-value.
	comp_max_interval(T, n, &ymin, &zmax);

	// Divide the interval into smaller intervals. Each smaller interval
	// contains one eigen-value.
	div_intervals(T, n, ymin, zmax, interval, &numInterval);

	// Compute the eigen-values one by one.
	for(i=0;i<numInterval;i++)
		eigenVal[i] = bisect_method(T, n, interval[i].y, interval[i].z);
	*numEigen = numInterval;
}


