/***********************************************************
 * This program demonstrate an example of root finding by using
 * Newton's method.
 */
#include <stdio.h>
#include <math.h>

#define EPSILON   0.000001


/*---------------------------------------------------------------------------
 * The function f(x, y) for the root-finding problem.
 *  f(x,y) = x^2 -5y -6
 */
double f(double x, double y)
{
    return (x*x-5.0*y - 6.0);
}

/*---------------------------------------------------------------------------
 * The function g(x, y) for the root-finding problem.
 *  g(x,y)= 3x - y^2  - 8
 */
double g(double x, double y)
{
    return (3.0*x-y*y - 8.0);
}

/*-------------------------------------------------
 * Derivatives of f() and g().
 */
double fx(double x, double y)
{
	return 2.0*x;
}

double fy(double x, double y)
{
	return (-5.0);
}

double gx(double x, double y)
{
	return (3.0);
}

double gy(double x, double y)
{
	return (-2.0*y);
}




/*--------------------------------------------------------------------------
 * The 2D Newton's method.
 *     *x, *y: the initial values of the root.
 */
void  newton2D(double *x, double *y)
{
  int        i=0;
  double  err, xnew, ynew, xold, yold;
  double  fdx, fdy, gdx, gdy, Delta;
  double  fold, gold;
  FILE   *fp;
  
  //Set up initial conditions.
  xold = *x;
  yold = *y;
//Function values & derivatives of the functions
  fold = f(xold, yold);
  gold = g(xold, yold);
  fdx = fx(xold, yold);
  fdy = fy(xold, yold);
  gdx = gx(xold, yold);
  gdy = gy(xold, yold);
//Determinant of the Jacobian matrix, no negative value.
  Delta = fdx*gdy - fdy*gdx;
//Update the root, Don't miss the minus operator!
  xnew = xold - (gdy*fold - fdy*gold)/Delta;
  ynew = yold - (-gdx*fold + fdx*gold)/Delta;
//Compute the 1-norm of the error vector.
  err = fabs(xnew - xold) + fabs(ynew-yold);
  fp = fopen("newton2D_result.txt","w");
  fprintf(fp, "      i          xn                   yn                    error\n");
  fprintf(fp,"------------------------------------------------------------\n");
  fprintf(fp,"     %3d \t%lf \t%lf \t---\n", i, xold, yold);
  i = 1;
  while(err>EPSILON){
   fprintf(fp,"     %3d \t%lf \t%lf \t%lf\n", i, xnew, ynew, err);

   xold = xnew; // Save current value.
   yold = ynew;
   //Update the root.
   fold = f(xold, yold);
   gold = g(xold, yold);
   fdx = fx(xold, yold);
   fdy = fy(xold, yold);
   gdx = gx(xold, yold);
   gdy = gy(xold, yold);

   Delta = fdx*gdy - fdy*gdx;

   xnew = xold - (gdy*fold - fdy*gold)/Delta;
   ynew = yold - (-gdx*fold + fdx*gold)/Delta;
   //Compute the 1-norm of the error vector.
   err = fabs(xnew - xold) + fabs(ynew-yold);
    i ++;
  }   
  fprintf(fp,"------------------------------------------------------------\n"); 
  fprintf(fp,"     %3d \t%lf \t%lf \t%lf\n", i, xnew, ynew, err);
 *x = xnew;
 *y = ynew;
}


/*----------------------------------------------------------------------------
  * The main procedure.
  */
int main(int argc, char **argv)
{
   double    x, y;


   fprintf(stderr," Input the 2 initial values x0, y0= ");
   fscanf(stdin,"%lf %lf", &x, &y);
   newton2D(&x, &y);
   fgetchar();
}
