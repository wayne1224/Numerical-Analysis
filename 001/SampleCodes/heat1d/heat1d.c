/*********************************************************
 * In this file an 1D heat transfer solver is presented.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//define num. of grid-points
#define MAXPNT   21
#define EPSILON  0.0001
#define a 0.0
#define b MAXPNT-1
#define h 1.0

double  t[MAXPNT];



/*-----------------------------------------------------
 * Set Bnd cond.
 */
void set_bound_cond(void)
{
	t[0] = 20.0;
	t[MAXPNT-1] = 100.0;
}


/*-------------------------------------------------------
 * Solve heat equation by using gauss-seidel method
 */
void gauss_seidel(void)
{
	int     i, k;
	double  err, y, infNorm;

	err = 1.0;
	for(i=0;i<MAXPNT;i++) t[i] = 0.0;
	k = 0;
	while(err>EPSILON){
		set_bound_cond();
		infNorm = 0.0;
		for(i=1;i<MAXPNT-1;i++){	  
          y = t[i];
		  t[i] = (t[i-1]+t[i+1])/(2*h*h);
		  if(fabs(t[i]-y)>infNorm)
			  infNorm = fabs(t[i]-y);
		}
		err = infNorm;
		k ++;
	}
}


void main(int argc, char **argv)
{
	gauss_seidel();
}


