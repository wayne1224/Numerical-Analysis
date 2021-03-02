#include <stdio.h>
#include <stdlib.h>
#include <math.h>



/*-----------------------------------------------------
 * Procedure to allocate space for an n-dimensional
 * vector. Return the pointer to the vector.
 */
double *alloc_vec(int n)
{
	double *t;

	t = (double *) malloc(sizeof(double)*n);
	return(t);
}


/*---------------------------------------------------
 * Procedure to initialize an n-dimensional vector.
 */
void init_vec(double *y, int n)
{
	int  i;

	/*--- init. vec. = 1. ---*/
	for(i=0;i<n;i++) y[i] = 1.0;
	y[0] = 0.0;
	y[0] = 0.8;
	y[1] = 0.3;
	y[2] = 0.4;
}


/*--------------------------------------------------
 * Compute inner product <a, b>, where a and b are 
 * n-dimensional vectors.
 */
double  inner_product(double *a, double *b, int n)
{
	int   i;
	double sum;

	sum = 0.0;
	for(i=0;i<n;i++)
		sum += a[i]*b[i];
	return (sum);
}


/*----------------------------------------------------------
 * Procedure to compute the outer product of two vector.
 *   A = x*y^T, where A is an NxN matrix, and x and y are
 * n dimensional vectors.
 */
void outer_product(double **A, double *x, double *y, int n)
{
	int   i, j;

	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			A[i][j] = x[i]*y[j];
}

/*-------------------------------------------------------------
 * Compute x= c*y, where x and y are N-dimensional vec. and c is 
 * a scalar.
 */
void scalar_vec_mult(double *x, double c, double *y, int n)
{
   int  i;

   for(i=0;i<n;i++) x[i] = c*y[i];
}

/*-------------------------------------------------------------
 * Normalize a vector. If the vector norm is 0, do nothing.
 */
void normalize_vec(double *a, int n)
{
	int    i;
	double sum;

	sum = 0.0;
	for(i=0;i<n;i++)
		sum += a[i]*a[i];
	if(sum==0.0) return;
	sum = sqrt(sum);
	for(i=0;i<n;i++)
		a[i] = a[i]/sum;
}

/*--------------------------------------------------------
 * Procedure to compute the 2-norm of a vector.
 */
double vec_norm(double *a, int n)
{
	double   sum;
	int      i;

	sum = 0.0;
	for(i=0;i<n;i++)
		sum += a[i]*a[i];
	sum = sqrt(sum);
	return(sum);
}

/*-----------------------------------------------------
 * Copy a vector from source to destination.
 *    *src: source 
 *    *dst: destination.
 */
void copy_vec(double *dst, double *src, int n)
{
	int   i;

	for(i=0;i<n;i++) dst[i] = src[i];
}


/*--------------------------------------------------------
 * Procedure to print out a vector.
 */
void print_vec(double *x, int n)
{
	int i;

    for(i=0;i<n;i++){
		fprintf(stderr," %lf", x[i]);
	}
	fprintf(stderr,"\n--------------------------\n");
}


/*#####################################################
  # Matrix-related computations.
  #####################################################
  */

/*--------------------------------------------------------
 * Procedure to allocate an n X n matrix and return the 
 * pointer to the matrix.
 */
double **alloc_mtx(int n)
{
	double  **B;
	int     i;

  B = (double **) malloc(sizeof(double *)*n);
  for(i=0;i<n;i++) 
	B[i] = (double *) malloc(sizeof(double)*n);
  return (B);
}




/*---------------------------------------------------------
 * Compute a = A*b, where a and b are vectors and A an n x n
 * matrix.
 */
void mtx_vec_mult(double *a, double **A, double *b, int n)
{
	int  i, j;
	
	for(i=0;i<n;i++){
		a[i] = 0.0;
		for(j=0;j<n;j++)
			a[i] += A[i][j]*b[j];
	}
}

/*---------------------------------------------------------
 * Compute a = A^t*b, where a and b are vectors and A an n x n
 * matrix.
 */
void transpose_mtx_vec_mult(double *a, double **A, double *b, int n)
{
  int   i, j;

  for(i=0;i<n;i++){
	  a[i] = 0;
	  for(j=0;j<n;j++) a[i] += A[j][i]*b[j];
  }
}



/*-------------------------------------------------------------
 * Compute B = c*A, where B and A are NxN matrices and c is 
 * a scalar.
 */
void scalar_mtx_mult(double **B, double c, double **A, int n)
{
	int i, j;

	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			B[i][j] = c*A[i][j];
}

/*--------------------------------------------------------
 * Procedure to compute matrix-matrix multiplication.
 *   A = B + C.
 */
void mtx_mtx_add(double **A, double **B, double **C, int n)
{
	int  i, j;

	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			A[i][j] = B[i][j] + C[i][j];
}


/*------------------------------------------------------------
 * Matrix-matrix multiplication.
 *   A = B*C.
 */
void mtx_mtx_mult(double **A, double **B, double **C, int n)
{
	int   i, j, k;

	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			A[i][j] = 0.0;
			for(k=0;k<n;k++)
				A[i][j] += B[i][k]*C[k][j];
		}
	}
}


/*------------------------------------------------------------
 * Procedure to make a matrix identity.
 */
void make_identity_mtx(double **A, int n)
{
	int i, j;

	for(i=0;i<n;i++){
		A[i][i] = 1.0;
		for(j=i+1;j<n;j++)
			A[i][j] = A[j][i] = 0.0;
	}
}


/*-------------------------------------------------------------
 * Procedure to compute residual vector  r[] = x[] - u*y[].
 *    u: the eigen value.
 */
void comp_residual(double *r, double *x, double u, double *y, int n)
{
	int   i;

	for(i=0;i<n;i++)
		r[i] = x[i] - u*y[i];
}

/*------------------------------------------------------
 * Procedure to print out a matrix.
 */
void print_mtx(double **A, int n)
{
	int i, j;

    for(i=0;i<n;i++){
		fprintf(stderr,"\n");
		for(j=0;j<n;j++)
			fprintf(stderr,"%lf ", A[i][j]);
	}
	fprintf(stderr,"\n--------------------------\n");
}


