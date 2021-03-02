/*************************************************************
* In this sample program. the Jocobi method is used to compute the eigenvalues of
* a symmetric matrix.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "definition.h"
#include "vec_mtx.h"

/*--- Declare the target matrix and its eiegn-vector matrix. ---*/
double   **A;
double   **P;

/*---  Declare eiegen value vector. ---*/
double   *eigenVec;

/*--- Declare test matrix for the orthogonality between eigen-vectors. ---*/
//orMtx[i][j] = inner_product(vi, vj); Should be a diagonal matrix.
double   **orMtx;

/*---  Declare  vector of residual norms. ---*/
double   *errVec;

/*------------------------------------------------------
 * Create a Hilbert matrix.
 *   A: the target mtx,
 *   n: dimension of matrix.
 */
void Hilbert_mtx(double **A, int n)
{
	int  i, j;

	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			A[i][j] = 1.0/(i+j+1.0);
}

/*-------------------------------------------------------------
 * Test matrix from Wiki
 */
void test_4_by_4_mtx(double **A)
{
    A[0][0]  =  4.0;
	A[1][0] = A[0][1] =  -30.0;
	A[1][1] = 300.0;
	A[2][0] = A[0][2] = 60.0;
	A[2][1] = A[1][2] = -675.0;
	A[2][2] = 1620.0;
	A[3][0] = A[0][3] = -35.0;
	A[3][1] = A[1][3] = 420.0;
	A[3][2] =  A[2][3] = -1050.0;
	A[3][3] = 700.0;
}

/*------------------------------------------------------
 * Create a symmetric  matrix.
 *   A: the target mtx,
 *   n: dimension of matrix.
 */
void symmetric_mtx(double **A, int n)
{
	int  i, j;


    for(i=0;i<n;i++)
		for(j=i;j<n;j++)
			A[i][j] = A[j][i] = (double)(rand()%50);
}


/*--------------------------------------------------------
 * Procedure to retrieve eigen values and eigen vectors
 * from a diagonal mtx and its eigen-vector matrix.
 */
void   retrieve_eigen_values(double **A, double **P, double *v, int n)
{
   int   i, j;

   for(i=0;i<n;i++) v[i] = A[i][i];
   // transpose the eigen-vector matrix.
   transpose_mtx(P, n);

   // Test the orthogonality of the eigen-vectors.
   orMtx = alloc_mtx(n);
   for(i=0;i<n;i++)
	   for(j=0;j<n;j++)
		   orMtx[i][j] = inner_product(P[i], P[j], n);

	fprintf(stderr,"Orthogonality between eigen-vectors:\n");
    for(i=0;i<n;i++){
		fprintf(stderr,"\n");
		for(j=0;j<n;j++)
			fprintf(stderr,"%lf ", orMtx[i][j]);
	}
	fprintf(stderr,"\n--------------------------\n");
}


/*---------------------------------------------------------
 * Procedure to compute the norms of
 *   ||A*x-u*x||, where A is the matrix, x the eigen vector,
 * and u the eigen value.
 * Compute this norm for all eigen values and their eigen vectors.
 * Keep them in err[].
 * Input:
 *    A: the mtx,
 *    P: mtx, each row = an eigen-vector.
 *    v[]: array of eigen-values.
 *    n: matrix dimension.
 */
void comp_err_vec(double **A, double **P, double *v, int n)
{
	int     i;
	double  *a, *b;

	a = alloc_vec(n);
	b = alloc_vec(n);
	errVec = alloc_vec(n);

	for(i=0;i<n;i++){
		//a = A*x.
      mtx_vec_mult(a, A, P[i], n);
	  // b = a - u*x
      comp_residual(b, a, v[i], P[i], n);
	  // err = ||b||
	  errVec[i] = vec_norm(b, n);
	}
	fprintf(stderr," Residual of A*x - u*x =\n");
	for(i=0;i<n;i++)
		fprintf(stderr," %lf", errVec[i]);
	fprintf(stderr,"\n-------------------------------------------------\n");
}


/*--------------------------------------------------------
 * Procedure to retrieve the eigen values and sort them 
 */

/*----------------------------------------------------
 * The main procedure
 */
int main(int argc, char **argv)
{
	int  i, j, n;

	n = 4;
	//Allocate mtx space.
    P = alloc_mtx(n);
    A = alloc_mtx(n);

	/*--- Initialize A[]. ---*/
	//Test case from Wiki
	n = 4;
    test_4_by_4_mtx(A);

  //  Hilbert_mtx(A, n);


  // 	symmetric_mtx(A, n);
   
	fprintf(stderr,"A[][]=\n");
    for(i=0;i<n;i++){
		fprintf(stderr,"\n");
		for(j=0;j<n;j++)
			fprintf(stderr,"%lf ", A[i][j]);
	}
	fprintf(stderr,"\n--------------------------\n");

	// Call Jacobian method to compute the eigen-values
	// and eigen-vectors.
    jacobian_method(A, P, n);
	fprintf(stderr,"After the diagonalization, A[][]=\n");
    for(i=0;i<n;i++){
		fprintf(stderr,"\n");
		for(j=0;j<n;j++)
			fprintf(stderr,"%lf ", A[i][j]);
	}
	fprintf(stderr,"\n--------------------------\n");

	fprintf(stderr,"The eigen-vector matrix P[][]=\n");
    for(i=0;i<n;i++){
		fprintf(stderr,"\n");
		for(j=0;j<n;j++)
			fprintf(stderr,"%lf ", P[i][j]);
	}
	fprintf(stderr,"\n--------------------------\n");

	/*---- Retrieve the eigen values. ----*/
	eigenVec = alloc_vec(n);
    retrieve_eigen_values(A, P, eigenVec, n);
  	fprintf(stderr,"EienValues= ");
    for(i=0;i<n;i++){
    	fprintf(stderr,"%lf ", eigenVec[i]);
	}
	fprintf(stderr,"\n--------------------------\n"); 
	fprintf(stderr,"The eigen-vectors are::\n");
    for(i=0;i<n;i++){
		fprintf(stderr,"x[%2d]=", i);
		for(j=0;j<n;j++)
			fprintf(stderr,"%lf ", P[i][j]);
	    fprintf(stderr,"\n--------------------------\n");

	}
	// Test the results
  	symmetric_mtx(A, n); //Re-generate the input matrix.
	 test_4_by_4_mtx(A);
    comp_err_vec(A, P, eigenVec, n);

	//Pause
	getchar();
}
