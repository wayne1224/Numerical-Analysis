#include <stdio.h>
#include <stdlib.h>
#include <math.h>



double   **A;
double   **L, **U;

/*------------------------------------------------------
 * Create a symmetric linear system.
 *   A: the coef. mtx,
 *   b: the right hand side.
 *   n: dimension of matrix=3
 */
void final_exam_mtx(double **A, double *b, int n)
{
   A[0][0] = 6.0; A[1][1] = 5.0; A[2][2] = 4.0;
   A[0][1] = A[1][0] = 3.0;
   A[0][2] = A[2][0] = 2.0;
   A[1][2] = A[2][1] = 1.0;
   b[0] = 11.0;
   b[1] = 9.0;
   b[2] = 7.0;
}

/*-------------------------------------------------------------------
 * Procedure to perform Doolittle LU decomposition.
 */
void doolittle(double **A, double **L, double **U, int n)
{
	double   temp;
    int   i, j, k;


	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			L[i][j] = U[i][j] = 0.0;

	for(i=0;i<n;i++){
		L[i][i] = 1.0;
	/*---- Compute the upper triangular matrix.
		  U_ij = A_ij - sum(L_ik*U_kj), k=0..i-1
		  */
		for(j=i;j<n;j++){
			temp = 0.0;
			for(k=0;k<=i-1;k++)
				temp = temp + L[i][k]*U[k][j];
			U[i][j] = A[i][j] - temp;
		}
		/*----Compute the lower triangular matrix. 
		      L_ij = A_ij - sum(L_ik*U_kj), k=0..i-1
			  */
		for(j=i+1;j<n;j++){
			temp = 0.0;
			for(k=0;k<=i-1;k++)
				temp = temp + L[j][k]*U[k][i];
			L[j][i] = (A[j][i]-temp)/U[i][i];
		}


	}
}

/*-----------------------------------------------------------
 * Forward substitution for solving lower triangular system.
 * L[][]: the lower triangular coefficient matrix,
 * y[]: the unknown
 * b[]: the right hand side.
 * n: diemension.
 */
void forward_substitute(double **L, double *y, double *b, int n)
{
	int  i, j;

	for(i=0;i<n;i++){
		y[i] = b[i]/L[i][i];
		for(j=i+1;j<n;j++)
			b[j] = b[j] - L[j][i]*y[i];
	}
}


/*-----------------------------------------------------------
 * Backward substitution for solving lower triangular system.
 * U[][]: the lower triangular coefficient matrix,
 * x[]: the unknown
 * y[]: the right hand side.
 * n: diemension.
 */
void backward_substitute(double **U, double *x, double *y, int n)
{
	int  i, j;

	for(i=n-1;i>=0;i--){
		x[i] = y[i]/U[i][i];
		for(j=i-1;j>=0;j--)
			y[j] = y[j] - U[j][i]*x[i];
	}
}


int main(int argc, char **argv)
{
	int  i, j;
	double  b[3], y[3], x[3];

    L = (double **) malloc(sizeof(double *)*3);
	for(i=0;i<3;i++) 
		L[i] = (double *) malloc(sizeof(double)*3);
    U = (double **) malloc(sizeof(double *)*3);
	for(i=0;i<3;i++) 
		U[i] = (double *) malloc(sizeof(double)*3);
    A = (double **) malloc(sizeof(double *)*3);
	for(i=0;i<3;i++) 
		A[i] = (double *) malloc(sizeof(double)*3);
/*
	A[0][0] = 6.0;
	A[0][1] = 3.0;
	A[0][2] = 2.0;
	A[1][0] = 3.0;
	A[1][1] = 2.0;
	A[1][2] = 1.5;
	A[2][0] = 2.0;
	A[2][1] = 1.5;
	A[2][2] = 1.2;
*/
  final_exam_mtx(A, b, 3);

	doolittle(A, L, U, 3);

	fprintf(stderr,"A[][]=\n");
    for(i=0;i<3;i++){
		fprintf(stderr,"\n");
		for(j=0;j<3;j++)
			fprintf(stderr,"%lf ", A[i][j]);
	}
	fprintf(stderr,"\n--------------------------\n");
	fprintf(stderr,"L[][]=\n");
    for(i=0;i<3;i++){
		fprintf(stderr,"\n");
		for(j=0;j<3;j++)
			fprintf(stderr,"%lf ", L[i][j]);
	}
	fprintf(stderr,"\n--------------------------\n");
	fprintf(stderr,"U[][] = \n");
    for(i=0;i<3;i++){
		fprintf(stderr,"\n");
		for(j=0;j<3;j++)
			fprintf(stderr,"%lf ", U[i][j]);
	}
	fprintf(stderr,"\n--------------------------\n");
//forward substitution for solving Ly=b
	forward_substitute(L,y,b, 3);
	fprintf(stderr,"y[]= ");
	for(i=0;i<3;i++) fprintf(stderr,"%lf ", y[i]);
	fprintf(stderr,"\n---------------------------------------\n");

	//backward  substitution for solving Ux=y
	backward_substitute(U, x, y, 3);
    fprintf(stderr,"x[]= ");
	for(i=0;i<3;i++) fprintf(stderr,"%lf ", x[i]);
	fprintf(stderr,"\n---------------------------------------\n");
	getchar();
}
