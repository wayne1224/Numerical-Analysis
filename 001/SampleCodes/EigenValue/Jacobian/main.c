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
	//Add shift matrix
//	for(i=0;i<n;i++) A[i][i] += 1.0/(2.0*n-1.0);
//	for(i=0;i<n;i++) A[i][i] += 1.0;
}

/*------------------------------------------------------
 * Create a symmetric  matrix.
 *   A: the target mtx,
 *   n: dimension of matrix.
 */
void symmetric_mtx(double **A, int n)
{
	int  i, j;

    srand(0);
    for(i=0;i<n;i++)
		for(j=i;j<n;j++)
			A[i][j] = A[j][i] = (double)(rand()%97);
	
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
 * Create a symmetric linear system.
 *   A: the coef. mtx,
 *   b: the right hand side.
 *   n: dimension of matrix=3
 */
void final_exam_mtx(double **A,  int n)
{
    A[0][0] = 6.0; A[1][1] = 5.0; A[2][2] = 4.0;
   A[0][1] = A[1][0] = 3.0;
   A[0][2] = A[2][0] = 2.0;
   A[1][2] = A[2][1] = 1.0;

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
	int     i, k;
	double  *a, *b, *t;

	a = alloc_vec(n);
	b = alloc_vec(n);
	errVec = alloc_vec(n);
	t = alloc_vec(n);

	for(i=0;i<n;i++){
		for(k=0;k<n;k++) t[k] = P[i][k];
		//a = A*x.
      mtx_vec_mult(a, A, t, n);
	  // b = a - u*x
      for(k=0;k<n;k++) b[k] = a[k] - v[i]*t[k];
	  // err = ||b||
	  errVec[i] = vec_norm(b, n);
	}
	fprintf(stderr," Residual of A*x - u*x =\n");
	for(i=0;i<n;i++)
		fprintf(stderr," %lf", errVec[i]);
	fprintf(stderr,"\n-------------------------------------------------\n");
}




/*----------------------------------------------------
 * The main procedure
 */
int main(int argc, char **argv)
{
	int  i, j, n, k;
	int  max, min;
	double condNum;

	n = 3;
	//Allocate mtx space.
    P = alloc_mtx(n);
    A = alloc_mtx(n);

	/*--- Initialize A[]. ---*/

//  Hilbert_mtx(A, n);
   final_exam_mtx(A, 3);
 //  symmetric_mtx(A, n);
 //    test_4_by_4_mtx(A);

	fprintf(stderr,"A[][]=\n");
    for(i=0;i<n;i++){
		fprintf(stderr,"\n");
		for(j=0;j<n;j++)
			fprintf(stderr,"%lf ", A[i][j]);
	}
	fprintf(stderr,"\n--------------------------\n");

	// Call Jacobian method to compute the eigen-values
	// and eigen-vectors.
    k = jacobian_method(A, P, n);
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
	fprintf(stderr,"No. of iterations= %d\n",k);
  	fprintf(stderr,"EienValues= ");
	max = min =0;
    for(i=0;i<n;i++){
    	fprintf(stderr,"%lf ", eigenVec[i]);
		if(fabs(eigenVec[i])>fabs(eigenVec[max])) max = i;
	    if(fabs(eigenVec[i])<fabs(eigenVec[min])) min = i;
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
  //	Hilbert_mtx(A, n); //Re-generate the input matrix.
//	symmetric_mtx(A, n);	
 //   test_4_by_4_mtx(A);
	final_exam_mtx(A, 3);
    comp_err_vec(A, P, eigenVec, n);
  //Compute the condition number
	condNum = fabs(eigenVec[max]/eigenVec[min]);
	fprintf(stderr, "max eigenvalue  = %.14f\n", eigenVec[max]);
	fprintf(stderr, "min eigenvalue = %.14f\n", eigenVec[min]);
	fprintf(stderr, "Condition number = %.14f\n", condNum);

	//Pause
	getchar();
}
