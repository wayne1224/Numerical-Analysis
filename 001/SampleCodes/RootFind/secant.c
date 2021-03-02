/***********************************************************
 * This program demonstrate an example of root finding by using
 * secant method.
 */
#include <stdio.h>
#include <math.h>

#define EPSILON   0.000001

double   a, b;


/*---------------------------------------------------------------------------
 * The function fo the root-finding problem.
 */
double f(double x)
{
    return (x*x-5.0*x - 6.0);
}




/*--------------------------------------------------------------------------
 * The bisection method.
 *     a, b: the initial values of the root.
 *  New value = 
 */
double  secant(double a, double b)
{
  int     i=1;
  double  err, x0, x1, xnew, f1, f0;
  FILE *fp;

  
  /*--- Set up initial conditions. ----*/
  x0= a;
  x1 = b;
  f1 = f(x1);
  f0 = f(x0);
  fp = fopen("secant_result.txt","w");
  // Print out table header.
  fprintf(fp, "     i            xn               error\n");
  fprintf(fp,"------------------------------------------------------------\n");

  //Start the root finding process.
  err = fabs(x1-x0);
  while(err>EPSILON){
   xnew = x1 - f1*(x1-x0)/(f1-f0);  // New value.
   fprintf(fp,"     %3d \t%lf \t%lf\n", i, xnew, err);
   // Update the two values and function values.
   x0 = x1;
   x1 = xnew;
   f0 = f1;
   f1 = f(x1);
   //Update error.
   err = fabs(x1-x0);
   i ++;
  }   
  fprintf(fp,"     %3d \t%lf \t%lf\n", i, xnew, err);

  fprintf(fp,"------------------------------------------------------------\n"); 
   fprintf(fp,"Root= %lf, error=%lf\n", xnew, f(xnew));

  return(xnew);
}


/*----------------------------------------------------------------------------
  * The main procedure.
  */
int main(int argc, char **argv)
{
   double    x0, x1;
   double    r;


   fprintf(stderr," Input the initial values x0, x1= ");
   fscanf(stdin,"%lf %lf", &x0, &x1);
   r = secant(x0, x1);
   fprintf(stderr,"Root= %lf, error=%lf\n", r, f(r));
   fgetchar();
}
