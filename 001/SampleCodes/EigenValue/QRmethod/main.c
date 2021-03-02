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
double   **A, **Q;
double   *eigenVals;
int      numEigen; 
int      n; // dimension of the system.

void Hilbert_mtx(double **, int);

/*------------------------------------------------------
 * Create a Hilbert linear system.
 *   A: the coef. mtx,
 *   b: the right hand side.
 *   n: dimension of matrix.
 */
void Hilbert_mtx(double **A, int n)
{
	int  i, j;

	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			A[i][j] = 1.0/(i+j+1.0);
		}
	}
}



/*------------------------------------------------------
 * Create a symmetric linear system.
 *   A: the coef. mtx,
 *   b: the right hand side.
 *   n: dimension of matrix.
 */
void symmetric_mtx(double **A, int n)
{
	int  i, j;

	for(i=0;i<n;i++){
		for(j=i;j<n;j++){
			A[i][j] = A[j][i] = rand()%10;
		}
	}
}

/*-----------------------------------------------------
 * Procedure to create a random matrix.
 */
void random_mtx(double **A, int n)
{
	int i, j;

	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			A[i][j] = rand() % 20;
}


/*----------------------------------------------------
 * The main procedure
 */
int main(int argc, char **argv)
{
    int i;

	n = 6;

    A = alloc_mtx(n);
    Q = alloc_mtx(n);
    make_identity_mtx(Q, n);
    eigenVals = (double *) malloc(sizeof(double)*n);

//	symmetric, or random;
//	random_mtx(A, n);
	symmetric_mtx(A, n);
//	Hilbert_mtx(A, n);
// Print out the initial linear system
	fprintf(stderr,"A[][]=\n");
    print_mtx(A, n);

// Print out the QR decomposition results.
	Hessenberg_form(A, Q, n);

	fprintf(stderr," Hessenberg form, A[][]=\n");
    print_mtx(A, n);

	fprintf(stderr," The orthogonal mtx P[][]=\n");
	print_mtx(Q, n);
    tri_diag_mtx_eigen(A, n, eigenVals, &numEigen);
	fprintf(stderr," The eigen values are:\n");
    for(i=0;i<n;i++) fprintf(stderr, "%14.8f ", eigenVals[i]);
	fprintf(stderr,"\n");
}
