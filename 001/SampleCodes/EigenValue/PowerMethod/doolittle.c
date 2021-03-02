#include <stdio.h>
#include <stdlib.h>
#include <math.h>




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
	/*--- Verify the solution.----*/
	for(i=0;i<n;i++)
		for(j=0;j<n;j++){
			temp = 0.0;
			for(k=0;k<n;k++)
				temp += L[i][k]*U[k][j];
			if(fabs(temp-A[i][j])>0.001)
				fprintf(stderr," Error in LU decomposition\n");
		}
}



