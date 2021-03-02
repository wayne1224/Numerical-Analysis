/*************************************************************
* In this sample program. the power method and inverse power
* method are used to compute the max and min eigenvalues of
* a matrix.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "definition.h"
#include "vec_mtx.h"

/*--- Declare the coef. matrix, the unkown vec. and the rhs.. ---*/
double   **A, *x, *b;
int      n; // dimension of the system.




/*------------------------------------------------------
 * Create a Hilbert linear system.
 *   A: the coef. mtx,
 *   b: the right hand side.
 *   n: dimension of matrix.
 */
void Hilbert_linear_system(double **A, double *b, int n)
{
	int  i, j;

	for(i=0;i<n;i++){
		b[i] = 0.0;
		for(j=0;j<n;j++){
			A[i][j] = 1.0/(i+j+1.0);
			b[i] += A[i][j];
		}
	}
}



/*------------------------------------------------------
 * Create a symmetric linear system.
 *   A: the coef. mtx,
 *   b: the right hand side.
 *   n: dimension of matrix.
 */
void symmetric_linear_system(double **A, double *b, int n)
{
	int  i, j;

	for(i=0;i<n;i++){
		for(j=i;j<n;j++){
			A[i][j] = A[j][i] = rand()%10;
		}
	}
	for(i=0;i<n;i++){
		b[i] = 0.0;
		for(j=0;j<n;j++) b[i] += A[i][j];
	}
}

/*------------------------------------------------------
 * Create a symmetric linear system.
 *   A: the coef. mtx,
 *   b: the right hand side.
 *   n: dimension of matrix=3
 */
void final_exam_mtx(double **A, double *b, int n)
{
   A[0][0] = A[1][1] = A[2][2] = 4.0;
   A[0][1] = A[1][0] = 1.0;
   A[0][2] = A[2][0] = 1.0;

   A[1][2] = A[2][1] = 2.0;
   b[0] = 4.0;
   b[1] = -1.0;
   b[2] = 3.0;
}

   
/*----------------------------------------------------
 * The main procedure
 */
int main(int argc, char **argv)
{

	n = 3;

    A = alloc_mtx(n);
    x = alloc_vec(n);
    b = alloc_vec(n);

//	Hilbert_linear_system(A, b, n);
//	symmetric_linear_system(A, b, n);
	final_exam_mtx(A, b, n);
// Print out the initial linear system
/*
	A[0][0] = 0; A[0][1] = -3; A[0][2] = 4;
	A[1][0] = 5; A[1][1] = 0; A[1][2] = -2;
    A[2][0] = 3; A[2][1] = 5; A[2][2] = 1;

	b[0] = -4; b[1] = 7; b[2] = 2;
*/
	fprintf(stderr,"A[][]=\n");
    print_mtx(A, n);
	fprintf(stderr,"b[]=\n");
	print_vec(b, n);
// Perform Gaussian elimination
    gauss_elm(A, b, n);
	fprintf(stderr," After forward elimination, A[][]=\n");
    print_mtx(A, n);
	//Perform backward substitution.
    back_substitute(A, x, b, n);
    //print out the results
	fprintf(stderr,"The solution x[]=\n");
	print_vec(x, n);
	getchar();
}
