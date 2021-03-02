/***************************************************************
 * This program solves the heat transfer problem of Hw4.
 * The governing equation: Laplace(T) = -s(x,y), s(0.5,0.5)=1.0, otherwise s(x,y)=0.0.
 * Domain: [0,1]x[0,1].
 *  Boundary conditions: T(0,y)=T(1,y)=T(x,0)=20.0. d(T)/d(n)=0 at y=1.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/*----- Define the error bound for stoping the program. ---*/
#define   EPSILON  0.0000001

/*-----The resolution of the regular grid, (N+1) by (N+1). -------------------------------*/
#define   N 20

//Define the extention of the domain
#define   WIDTH   1.0
#define   HEIGHT  1.0

//Boundary value of the Dirichlet BC
#define   BC_VALUE    20.0
//Derivative value of the Neumann BC
#define   BC_DERIVATIVE   -50

//Define the source value at the source
#define   SRC_VALUE    -10000.0

//The grid-points
double   T[N+1][N+1];

/*----------------------------------------------------------------------------------
 * Assign the Dirichlet BC at the  bottom of the domain.
 */
void  set_Dirichlet_BC(void)
{
	int   i, j;

	for(i=0;i<=N;i++) T[i][0] = BC_VALUE;
}

/*-----------------------------------------------------------------------------------
 * Procedure to set Neumann BC at the top boundary, j=N.
 *    d(T)/d(n) = b,  T[i][N] = b + h*T[i][N-1].
 */
void set_Neumann_BC(double w)
{
	int            i;
	double     h;

	h = WIDTH/N;
	for(i=0;i<=N;i++)  T[i][N] = T[i][N] + w*(T[i][N-1] + h*BC_DERIVATIVE-T[i][N]);
}

/*-------------------------------------------------------------------------------------
 * The source function -s(x,y).
 *      -s(WIDTH/2,HEIGHT/2) = -1.0;
 *     otherwise, s(x, y) = 0;
 *  Input:
 *      i, j: indices of the grid-point.
 *  Return the value of source.
 */
double   source(int i, int j)
{
	 if(i==N/2 && j==N/2) return (SRC_VALUE);
	 else return(0.0);
}


/*--------------------------------------------------------------------------------------
 * The SOR solver.
 *    1. Sweeping starts from T[1][1] to T[N-1][N-1].
 *    2. Neumann BC is applied upon T[i][N], the top boundary.
 *    3. Dirichlet BC is set upon the left, right and bottom boundaries.
 * Input:
 *    w: the relaxation weight.
 * Return:
 *    The infinite norm of T(k+1)-T(k).
 */
double SOR_iteration(double w)
{
    double   s, h, error, oldT, temp;
	int          i, j;

	/*-- One sweep.  ---*/
	error = 0.0;
	h = WIDTH/N;
	//The internal nodes
	for(i=1;i<=N-1;i++){
		for(j=1;j<=N-1;j++){
			oldT = T[i][j];
			s = source(i, j);
			T[i][j] = T[i][j] + (w/4.0)*(-s*h*h + (T[i-1][j]+T[i][j-1]+T[i][j+1]+T[i+1][j])-4.0*T[i][j]);
            temp =fabs(T[i][j]-oldT);
			if(temp>error) error = temp;
		}
	}
	//The left and the right boundaries
	i = 0;
	for(j=1;j<=N-1;j++){
		T[i][j] = BC_VALUE;
	}
	i = N;
	for(j=1;j<=N-1;j++){
		T[i][j] = BC_VALUE;
	}

	/*--- apply the Neumann BC. ---*/
	set_Neumann_BC(w);
	//Return the error, the infinite norm
	return(error);
}


/*----------------------------------------------------------------------------------------
 * The main procedure.
 */
int main(int argc, char **argv)
{
	int                 k, i, j;
	double          error, w;
	char              filename[64]="solution.txt";
	FILE            *fp;

	//Initialize the function values
    for(i=0;i<=N;i++)
		for(j=0;j<=N;j++)
			T[i][j] = 0.0;

	//Initialize the BC values.
    set_Dirichlet_BC();

	//Perform SOR iteration until the error is bounded.
	k=0;
	error = 999.0;
	w = 1.2;
	while(error>EPSILON){
		error = SOR_iteration(w);
		fprintf(stderr," k=%d, error=%lf\n", k, error);
		k++;
	}
	//Output the results
	fp = fopen(filename,"w");
	fprintf(fp,"%d %d %lf\n", N, N, WIDTH/N);
	for(i=0;i<=N;i++){
		for(j=0;j<=N;j++) fprintf(fp,"%lf ",T[i][j]);
        fprintf(fp,"\n");
	}
	fclose(fp);
	getchar();
}