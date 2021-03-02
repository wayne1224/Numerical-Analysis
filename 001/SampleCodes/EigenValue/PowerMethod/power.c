#include <stdio.h>
#include <math.h>

#include "definition.h"
#include "vec_mtx.h"



//L & U matrix for the matrix A[][]
double   **L, **U;
double   *y, *r, *b;


/*----------------------------------------------------------------
 * Procedure to solve a lower triangular system
 *    L*x = b.
 */
void forward_substitue(double **L, double *x, double *b, int n)
{
	int   i,  j;

	for(i=0;i<n;i++){
		x[i] = b[i]/L[i][i];
		for(j=i+1;j<n;j++)
			b[j] = b[j] - L[j][i]*x[i];
	}
}

/*----------------------------------------------------------------
 * Procedure to solve an upper triangular system
 *    U*x = b.
 */
void back_substitue(double **U, double *x, double *b, int n)
{
	int   i,  j;

	for(i=n-1;i>=0;i--){
		x[i] = b[i]/U[i][i];
		for(j=i-1;j>=0;j--)
			b[j] = b[j] - U[j][i]*x[i];
	}
}	   

/*------------------------------------------------------------
 * Procedure to solve L*U*x = y, L and U are lower & upper triangular
 * mtx, y is thr rhs, and x is the unknown.
 */
void solve_L_U_sys(double **L, double **U, double *x, double *y, int n)
{


	//Solve L*b = y for b[].
	forward_substitue(L, b, y, n);

	//Solve U*x = b for x[].
	back_substitue(U, x, b, n);


}


/*-------------------------------------------------------------
 * Power method to compute the max eigenvalue and eigenvector.
 *   A: an n X n mtx,
 *   x: the max eigen vector,
 *   val: the pointer to the max eigen value,
 *   n: dimension of the matrix.
 */
void power_method(double  **A, double *x, double *val, int n)
{
	int     k, i;
	double  err, u;


	y = alloc_vec(n);
	r = alloc_vec(n);
	init_vec(y, n);

	err = 1.0;
	k = 0;
	while(err>EPSILON && k<MAXSTEP){
		normalize_vec(y, n);

		fprintf(stderr,"y[]=");
		for(i=0;i<n;i++)
			fprintf(stderr," %lf ", y[i]);
		fprintf(stderr,"\n");

		mtx_vec_mult(x, A, y, n);

		u = inner_product(y, x, n);
		comp_residual(r, x, u, y, n);
		err = vec_norm(r, n);
		copy_vec(y, x, n);
		k ++;

		fprintf(stderr," k=%d, residual = %lf, eigenVal=%lf\n", k, err, u);
		fprintf(stderr,"x[]=");
		for(i=0;i<n;i++)
			fprintf(stderr," %lf ", x[i]);
		fprintf(stderr,"\n");

	}
	*val = u;
	normalize_vec(x, n);
}


/*-------------------------------------------------------------
 * Inverse Power method to compute the min eigenvalue and eigenvector.
 *   A: an n X n mtx,
 *   x: the min eigen vector,
 *   val: the pointer to the min eigen value,
 *   n: dimension of the matrix.
 */
void inverse_power_method(double  **A, double *x, double *val, int n)
{
	int     k, i;
	double  err, u;


	//Allocate L, U mtx space.
    L = alloc_mtx(n);
    U = alloc_mtx(n);
  // LU-decompose A[].
	doolittle(A, L, U, n);
	// Allocate temporary vector for backward substitution
	b = alloc_vec(n);

	init_vec(y, n);
	err = 1.0;
	k = 0;
	while(err>EPSILON && k<MAXSTEP){
		// Solve L*U*x = y to obtain new vector x[].
		normalize_vec(y, n);
		//Print y[]
		fprintf(stderr,"y[]=");
		for(i=0;i<n;i++)
			fprintf(stderr," %lf ", y[i]);
		fprintf(stderr,"\n");
		//Compute x = A^(-1)y
		solve_L_U_sys(L, U ,x , y, n); // x=A^(-1)y

	// Verify the solution
	mtx_vec_mult(b, A, x, n);
	comp_residual(r, b, 1.0, y, n);
	err = vec_norm(r, n);
	fprintf(stderr," ERR=%lf\n", err);

		u = inner_product(y, x, n);
		comp_residual(r, x, u, y, n);
		err = vec_norm(r, n);
		copy_vec(y, x, n);
		k ++;

		// print out intermediate results.
		fprintf(stderr," k=%d, residual = %lf, eigenvalue=%lf \n", k, err, u);
		fprintf(stderr,"x[]=");
		for(i=0;i<n;i++)
			fprintf(stderr," %lf ", x[i]);
		fprintf(stderr,"\n");
		
	}
	normalize_vec(y, n);
	solve_L_U_sys(L, U ,x , y, n); // x=A^(-1)y
	u = inner_product(y, x, n);

	normalize_vec(x, n);
	*val = 1.0/u;
}