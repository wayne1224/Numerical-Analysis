/***********************************************************
 * This program demonstrate an example of root finding by using
 * Newton's method.
 */
#include <stdio.h>
#include <math.h>

#define EPSILON   0.000000001

double   a, b;


/*---------------------------------------------------------------------------
 * The function fo the root-finding problem.
 * f(x) = x^2 - 6.0
 */
double f(double x)
{
    //return (x*x - 6.0);
    //Test case for double roots
	return ((x-2.0)*(x-2.0)*(x*x+x+1.0));
}


/*--------------------------------------------------------------------------
 * derivative of the target function.
 * fx(x) = 2*x;
 */
double fx(double x)
{
//  return(2.0*x);
   //Test case for double roots
	return (2.0*(x-2.0)*(2.0*x+1.0));
}

/*--------------------------------------------------------------------------
 * The bisection method.
 *     a: the initial value of the root.
 */
double  newton(double a)
{
  int     i=1;
  double  err, xnew, xold;
  FILE   *fp;

  fp = fopen("newton_result.txt","w");
  //Set up initial conditions.
  xold = a;
 // xnew = xold - f(xold)/fx(xold);
  //For double-root case
  xnew = xold - 1.414*f(xold)/fx(xold);
  err = fabs(xnew - xold);

  fprintf(fp, "     i            xn               error\n");
  fprintf(fp,"------------------------------------------------------------\n");
  fprintf(fp,"     %3d \t%15.10f \t ---\n", 0, xold);

  while(err>EPSILON){
   fprintf(fp,"     %3d \t%15.10f \t%15.10f\n", i, xnew, err);

   xold = xnew; // Save current value.
 // Compute new values.
 //  xnew = xold - f(xold)/fx(xold); 
   //For double-root case
   xnew = xold - 1.414*f(xold)/fx(xold);
   err = fabs(xnew - xold); // Compute the difference.
   i ++;
  }   
  fprintf(fp,"     %3d \t%15.10f \t%15.10f\n", i, xnew, err);

  fprintf(fp,"------------------------------------------------------------\n"); 
  fprintf(fp,"Root= %15.10f, f(x)=%15.10f, err=%15.10f\n", xnew, f(xnew), err);

  return(xnew);
}


/*----------------------------------------------------------------------------
  * The main procedure.
  */
int main(int argc, char **argv)
{
   double    r;


   fprintf(stderr," Input the initial value x0= ");
   fscanf(stdin,"%lf", &r);
   r = newton(r);
   fprintf(stderr,"Root= %lf, error=%lf\n", r, f(r));
   fgetchar();
}
