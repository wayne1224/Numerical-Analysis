/*****************************************************
 * Thsi file contains procedure to tranform a nxn matrix
 * into an upper Hessenberg matrix by using Householder
 * matrices.  
 *     
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "definition.h"
#include "vec_mtx.h"


// Temperary matrix & vector for similarity transformation.
double  **T, *x;

/*------------------------------------------------------
 * Create a Householder vector for a matrix. This vector
 * is used to create a Householder matrix H such that 
 * H*A makes the entries below A[i][j] of A[][j] be 0.
 * In solving linear system i==j. However, for producing
 * Hessenberg form, i=j+1.
 */
void create_H_vec(double **A, int i, int j, double *v, int n)
{
	int      k;
    double   norm2;

	for(k=0;k<i;k++) v[k] = 0.0;
	norm2 = 0.0;
	for(k=i;k<n;k++) norm2 += A[k][j]*A[k][j];
	norm2 = sqrt(norm2);

	if(A[i][j]>=0.0) v[i] = A[i][j] + norm2;
	else v[i] = A[i][j] - norm2;

	for(k=i+1;k<n;k++) v[k] = A[k][j];
}



/*--------------------------------------------------------
 * Procedure to compute A = A*P = A(I-2*V*V^t/(V^t*V))
 *  A= A + b*(A*V)*V^t = A + b*outer(A*V, V).
 */
void post_mult_P_mtx(double **A, double *v, int n)
{
	double   b;

	//b = -2.0/(V^t*V);
	b = -2.0/inner_product(v, v, n);
	// x = A*V;
	mtx_vec_mult(x, A, v, n);
	// T = x*V^t;
	outer_product(T, x, v, n);
	// T = b*T;
    scalar_mtx_mult(T, b, T, n);
	// A = A + T;
	mtx_mtx_add(A, A, T, n);
}


/*--------------------------------------------------------
 * Procedure to compute A = P*A = (I-2*V*V^t/(V^t*V))A
 *  A= A + b*V*(V^t*A) = A + b*V*(A^t*V)^t
 *   = A + b*outer(V, A^t*V).
 */
void pre_mult_P_mtx(double **A, double *v, int n)
{
	double   b;

	//b = -2.0/(V^t*V);
	b = -2.0/inner_product(v, v, n);
	// x = A^t*V;
	transpose_mtx_vec_mult(x, A, v, n);
	// T = V*x^t;
	outer_product(T, v, x, n);
	// T = b*T;
    scalar_mtx_mult(T, b, T, n);
	// A = A + T;
	mtx_mtx_add(A, A, T, n);
}


/*------------------------------------------------------
 * Procedure to do similarity transformation:
 *   A = PAP, 
 *   Q = Q*P.
 * Such that A becomes a Hessenberg matrix.
 */
void Hessenberg_form(double **A, double **Q, int n)
{
	int     j;
	double  *v;

	// Allocate temporary spaces.
	v = alloc_vec(n);
	T = alloc_mtx(n);
    x = alloc_vec(n);
	//Eliminate each column to make A[][] into upper Hessenberg
	//form.
	for(j=0;j<n-2;j++){
		// Create a vector v[] = A.j + c*e1, 
		// c=sign(A[j+1][j])*norm2(A[j+1][j]~A[n-1][j]).
		create_H_vec(A, j+1, j, v, n);
        
		// Compute A = P*A*P,
		post_mult_P_mtx(A, v, n);
        pre_mult_P_mtx(A, v, n);

		// Compute Q = Q*P
		post_mult_P_mtx(Q, v, n);
	}
}


 