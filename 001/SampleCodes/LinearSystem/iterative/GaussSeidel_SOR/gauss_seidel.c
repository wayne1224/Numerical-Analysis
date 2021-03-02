#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define  EPSILON    0.000001
#define  MAXN       3

double   **A, *b, *x, *y;
char     filename[32] = "output.dat";
FILE     *fp;


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

/*---------------------------------------------------------------
 * Procedure to compute the difference between two vectors.
 *   x, y: the two vectors
 *   N: dimension of the vectors.
 * Return the norm of the difference.
 */
double vec_diff(double *x, double *y, int N)
{
	double     sum;
	int        i;

	sum = 0.0;
	for(i=0;i<N;i++)
		sum = sum + (x[i]-y[i])*(x[i]-y[i]);
	sum = sqrt(sum);
	return (sum);
}


/*-------------------------------------------------------------------
 * Procedure to solve Ax=b by using Jacobian Method
 *   A is an n X n matrix.
 */
void jacobian(double **A, double *x, double *b, int n)
{
	double  err, sum;
    int     i, j, k;


	err = 1.0;
    
	//Guess an initial solution.
	for(i=0;i<n;i++) y[i] = x[i] = 0.0;
     
	k = 0; // number of iterations
	while(err>EPSILON){
		// Compute X^(n+1)[i]
		for(i=0;i<n;i++){
			sum = 0.0;
			for(j=0;j<n;j++)
				if(i!=j) sum = sum + A[i][j]*y[j];
			x[i] = (b[i]-sum)/A[i][i];
		}
		err = vec_diff(x, y, n);
		// Record new X[] in y[].
		for(i=0;i<n;i++) y[i] = x[i];
		k ++; // Increase num. of iterations
		fprintf(fp," k=%d,  err= %lf\n", k, err);
  	   fprintf(stderr," k=%d,  err= %lf\n x[]= ", k, err);
	   for(i=0;i<n;i++) fprintf(stderr," %lf ", x[i]);
	   fprintf(stderr,"\n");
	}	
}


/*-------------------------------------------------------------------
 * Procedure to solve Ax=b by using Gauss_Seidel Method
 *   A is an n X n matrix.
 */
void gauss_seidel(double **A, double *x, double *b, int n)
{
	double  err, sum;
    int     i, j, k;


	err = 1.0;
    
	//Guess an initial solution.
	for(i=0;i<n;i++) y[i] = x[i] = 0.0;
     
	k = 0; // number of iterations
	while(err>EPSILON){
		// Compute X^(n+1)[i]
		for(i=0;i<n;i++){
			sum = 0.0;
			for(j=0;j<n;j++)
				if(i!=j) sum = sum + A[i][j]*x[j];
			x[i] = (b[i]-sum)/A[i][i];
		}
		err = vec_diff(x, y, n);
		// Record new X[] in y[].
		for(i=0;i<n;i++) y[i] = x[i];
		k ++; // Increase num. of iterations
		fprintf(fp," k=%d,  err= %lf\n", k, err);
    	fprintf(stderr," k=%d,  err= %lf\n x[]= ", k, err);
		for(i=0;i<n;i++) fprintf(stderr," %lf ", x[i]);
		fprintf(stderr,"\n");
	}	
}


/*-------------------------------------------------------------------
 * Procedure to solve Ax=b by using SOR Method
 *   A is an n X n matrix.
 */
void SOR(double **A, double *x, double *b, int n, double w)
{
	double  err, sum;
    int     i, j, k;



	err = 1.0;
    
	//Guess an initial solution.
	for(i=0;i<n;i++) y[i] = x[i] = 0.0;
    x[1] = y[1] = 1.0;
	k = 0; // number of iterations
	while(err>EPSILON){
		// Compute X^(n+1)[i]
		for(i=0;i<n;i++){
			sum = 0.0;
			for(j=0;j<n;j++)
				if(i!=j) sum = sum + A[i][j]*x[j];
			x[i] = (1.0-w)*x[i]+ w*(b[i]-sum)/A[i][i];
		}
		err = vec_diff(x, y, n);
		// Record new X[] in y[].
		for(i=0;i<n;i++) y[i] = x[i];
		k ++; // Increase num. of iterations
		fprintf(fp," k=%d,  err= %lf\n", k, err);
    	fprintf(stderr," k=%d,  err= %lf\n x[]= ", k, err);
		for(i=0;i<n;i++) fprintf(stderr," %lf ", x[i]);
		fprintf(stderr,"\n");
	}	
}



int main(int argc, char **argv)
{
	int  i, j;

	// Allocate memory for coef matrx, rhs, and unkowns.
    x = (double *) malloc(sizeof(double)*MAXN);
	b = (double *) malloc(sizeof(double)*MAXN);
    A = (double **) malloc(sizeof(double *)*MAXN);
	y = (double *) malloc(sizeof(double));
	for(i=0;i<MAXN;i++) 
		A[i] = (double *) malloc(sizeof(double)*MAXN);

    // Use Hilter mtx as coef. mtx, and x[] = 1
	for(i=0;i<MAXN;i++){
		A[i][i] = 0.0;
		for(j=0;j<MAXN;j++){
		  if(i!=j) A[i][j] = 1.0/(i+j+1.0);
		}
		for(j=0;j<MAXN;j++)  // Diagonal dominant mtx
			 A[i][i] += A[i][j];
	}
	for(i=0;i<MAXN;i++){
		b[i] = 0.0;
		for(j=0;j<MAXN;j++)
			b[i] = b[i] + A[i][j];
	}
	final_exam_mtx(A, b, 3);
    // Open log file for recording
	fp = fopen(filename,"w");
    // Solve the system by using Jacobian method
 //   jacobian(A, x, b, MAXN); 
	// Solve the system by using Gauss-Seidel method
	gauss_seidel(A, x, b, MAXN);
//	SOR(A, x, b, MAXN, 1.20);
	fclose(fp);

	fprintf(stderr,"\n------------------------------------------------------------------\n");
	for(i=0;i<MAXN;i++)
		fprintf(stderr," x[%d]= %lf\n", i, x[i]);
	fprintf(stderr,"\n");
	fgetchar();
}
