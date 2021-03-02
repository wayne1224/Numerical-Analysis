/***********************************************************
 * This program demonstrate an example of root finding by using
 * bisection method.
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
    return (x*x-3.0*x+1.0);
}

/*--------------------------------------------------------------------------
 * The bisection method.
 *     a, b: the two ends of the interval. a<b.
 */
double  bisect(double a, double b)
{
  int i=0;
  double  errA, errB, errC, c;
  FILE    *fp;

  fp = fopen("bisect_result.txt", "w");
  if(f(a)>0.0 && f(b)>0.0){
    fprintf(stderr,"Bad initial interval. %lf - %lf.  Return.\n", a, b);
    return (a);
 }

   fprintf(fp, "     i      a            b             c                errA               errB             errC\n");
   fprintf(fp,"------------------------------------------------------------\n");
  while(fabs(b-a)>EPSILON){
       c = (a+b)/2.0;
       errA = f(a); errB = f(b); errC=f(c);
       i ++;
       fprintf(fp,"%3d %lf    %lf   %lf     %lf    %lf    %lf\n",
                               i, a, b, c, errA, errB, errC);
       if(errC==0.0) return (c);
       if(errA<=0.0 &&errC<=0.0||errA>=0.0&&errC>=0.0) a = c;
       else b = c;
   }
   fprintf(fp,"------------------------------------------------------------\n"); 
   fprintf(fp,"Root= %lf, error=%lf\n", c, f(c));
   fclose(fp);
   return(c);
}


/*----------------------------------------------------------------------------
  * The main procedure.
  */
int main(int argc, char **argv)
{
   double    r;
   char        c;


   fprintf(stderr," Input the initial interval a and b= ");
   fscanf(stdin,"%lf %lf", &a, &b);
   r = bisect(a, b);
 
   c = fgetc(stdin);
}
