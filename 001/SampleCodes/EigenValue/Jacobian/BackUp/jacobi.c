#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "definition.h"
#include "vec_mtx.h"

/*------------------------------------------------------
 * Procedure to print out a matrix.
 */
void print_mtx(double **A, int n)
{
	int i, j;

    for(i=0;i<n;i++){
		fprintf(stderr,"\n");
		for(j=0;j<n;j++)
			fprintf(stderr,"%12.8f ", A[i][j]);
	}
	fprintf(stderr,"\n--------------------------\n");
}


/*-------------------------------------------------------------
 * Jacobian method to compute the eigenvalues and eigenvectors.
 *   A: the input matrix,an n X n mtx,
 *   P: matrix of the eigen vectors,
 *   n: dimension of the matrix.
 * Matrix A will be diagonalized. The eigen values are at the
 * main diagonal.
 */
void jacobian_method(double  **A, double **P,  int n)
{
	double   **R;
	int      k, p, q;
	double   A_pq, c, s;
	double   theta, a, b;

	// Make P = I.
	make_identity_mtx(P, n);
	// Create a rotation matrix
	R = alloc_mtx(n);
	// Find the max off-diagonal entry.
	A_pq = max_off_diag_entry(A, &p, &q, n);
	k = 0;
	while(fabs(A_pq)>EPSILON && k<MAXSTEP){
		a = 2*A_pq;
		b = A[p][p]-A[q][q];
		if(b==0.0){
		  if(a>0.0) theta = PI/2.0;
		  else theta = 3.0*PI/2.0;
		}else{
		  theta = atan(a/b);
		}
		theta = theta/2.0;
		c = cos(theta);
		s = sin(theta);
        make_rotate_mtx(R, p, q, c, s, n);
        //P = P*R.
		post_mtx_mult(P, R, n); 
		// A = R'*A*R.
		post_mtx_mult(A, R, n);
		transpose_mtx(R, n);
		pre_mtx_mult(R, A, n);
		A[p][q] = A[q][p] = 0.0; // Enforce the max entry=0
    
		// Find the max off-diagonal entry for next iteration.
	    A_pq = max_off_diag_entry(A, &p, &q, n);
		k++;
    
		fprintf(stderr,"k=%d, matrix A[][]=\n", k);
		print_mtx(A, n);
	}

}

