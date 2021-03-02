/*************************************************************
* In this sample program, the Conjugate Grdaient Method is used to solve a linear
*  system Ax = b, where A is Symmetric Positive Definite (SPD).
*   S. K. Ueng, 2017/12/09
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "conjugate_gradient.h"
#include "vec_mtx.h"

/*--- Declare the coef. matrix, the unkown vec. and the rhs.. ---*/
double   **A, *x, *b;
int          n; // dimension of the system.


/*-------------------------------------------------------------------
 * Procedure to solve Ax=b by using SOR Method
 *   A is an n X n matrix.
 */
void SOR(double **A, double *x, double *b, int n, double w)
{
	double  err, sum, temp;
    int     i, j, k;
	
	err = 1.0;
  	//Guess an initial solution.
//	init_vec(x, n);
	    
	k = 0; // number of iterations
	while(err>EPSILON||k<n*2){
		// Compute X^(n+1)[i]
		err = 0.0;
		for(i=0;i<n;i++){
			temp = x[i];
			sum = 0.0;
			for(j=0;j<n;j++)
				 sum = sum + A[i][j]*x[j];
			//Update x[i]
			x[i] = x[i]+ w*(b[i]-sum)/A[i][i];
			temp = temp - x[i]; //Record the difference
			err += temp*temp; //Accumulate the norm
		}
	    err = sqrt(err);   //Compute the norm.
		k ++; // Increase num. of iterations
		fprintf(stderr," k=%d,  err= %lf\n", k, err);
    	fprintf(stderr," x[]= ");
		for(i=0;i<n;i++) fprintf(stderr," %lf ", x[i]);
		fprintf(stderr,"\n");
	}	
}

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
	int  t;

   srand(0);
	for(i=0;i<n;i++){
		for(j=i;j<n;j++){
			t = rand() % 10;
//			A[i][j] = A[j][i] = (double) (t+1.0);
			A[i][j] = A[j][i] = (i+j+1.0);
		}
 A[i][i] += 1.0;
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
    A[0][0] = 6.0; A[1][1] = 5.0; A[2][2] = 4.0;
   A[0][1] = A[1][0] = 3.0;
   A[0][2] = A[2][0] = 2.0;
   A[1][2] = A[2][1] = 1.0;
   b[0] = 11.0;
   b[1] = 9.0;
   b[2] = 7.0;
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

  //Hilbert_linear_system(A, b, n);
 //  symmetric_linear_system(A, b, n);
	final_exam_mtx(A, b, 3);
// Print out the initial linear system
   fprintf(stderr,"A[][]=\n");
    print_mtx(A, n);
    fprintf(stderr,"b[]=\n");
    print_vec(b, n);
// Print out the  results.
    conjugate_gradient(A, x, b, n);
    fprintf(stderr,"The solution x[]=\n");
    print_vec(x, n);

	//Try SOR
	/*
	fprintf(stderr,"\n *****************************SOR *******************\n");
    Hilbert_linear_system(A, b, n);
    SOR(A, x, b, n, 0.5);
	*/
    getchar();
}
