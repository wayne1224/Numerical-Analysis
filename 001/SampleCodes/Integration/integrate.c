/****************************************************************************
 * This sample program reveals the accuracy of integration algorithms.
 * The algorithms: Riemann sum, Trapezoid rules, Simpson rule.
 */
#include <stdio.h>
#include <math.h>

#define   N  80

//Sample points
double  f[N+1];
//integral range.
double  a=0.0, b=1.0;


/*------------------------------------------------------------------------------------------------------
 * The function to be integrated: f(x) = 1.0/(1+x^2).
 * The indefinite integral I(f(x)=F(x) =  arctan(x)
 */
double  my_func(double x)
{

   return(1.0/(1.0+x*x));
}


/*-----------------------------------------------------------------------------------------------------
 * The inetgration of f(x) = 1.0/(1+x^2).
 *  F(x) =arctan(x)
 */
double my_func_integral(double  a, double b)
{
   double     t;
   double     Fa, Fb;

    Fa = atan(a);
     Fb =atan(b);
   return(Fb-Fa);
}
   
/*------------------------------------------------------------------------------------------------------
 * Sample point generation.
 *    n: number of intervals.
 * Number of smplaes = n+1.
 */
void data_gen(int n)
{
   int         i;
   double   h;

   h = (b-a)/n;
   for(i=0;i<=n;i++)
      f[i] = my_func(i*h);
}

/*-----------------------------------------------------------------------------------------------------
 * Integration by Riemann sum
 */
double  Riemann(double f[], double h, int n)
{
   double  sum=0.0;
   int        i;

   for(i=0;i<n;i++) sum += f[i];
   return(sum*h);
}

/*--------------------------------------------------------------------------------------------------------
 * Integration by Trapezoid rule
 */
double Trapezoid(double f[], double h, int n)
{
   double    sum=0.0;
   int          i;

   sum =  f[0] + f[n];
   for(i=1;i<=n-1;i++)
      sum += 2.0*f[i];
   return(sum*h/2.0);
}

/*---------------------------------------------------------------------------------------------------------
 * Inetgration by using Simpson's rule.
 */
double   Simpson(double f[], double h, int n)
{
   int     i;
   double   sum;

   sum = f[0] + f[n];
   for(i=1;i<=n-1;i++){
      if(i%2==0.0) sum += 2.0*f[i];
      else sum += 4.0*f[i];
  }
   return(sum*h/3.0);
 }

/*----------------------------------------------------------------------------------------------------------
 * The main procedure 
 */
void main(int argc, char **argv)
{
  char      filename[]="integration.txt";
  int        numInterval;
  double   riemannSum, trapezoidSum, simpsonSum;
  double   exactSum, h;
  FILE     *fp;


  //Open the output file
  fp = fopen(filename, "w");
  //Compute the exact sum.
  exactSum = my_func_integral(a, b);
  fprintf(fp,"Exact integration = %lf\n", exactSum);
  fprintf(fp,"-----------------------------------------------------------------------------\n");
  //Try n= 10, 20, 30
  for(numInterval = 10;numInterval<=N;numInterval=numInterval+10){
     data_gen(numInterval);
    //Compute the numerical integrations
     h = (b-a)/numInterval;
     riemannSum = Riemann(f, h, numInterval);
     trapezoidSum = Trapezoid(f, h, numInterval);
     simpsonSum = Simpson(f, h, numInterval);
     fprintf(fp,"N= %d, Riemann= %lf, Trapezoid=%lf, Simpson=%lf\n", numInterval, riemannSum, 
                  trapezoidSum, simpsonSum);
 }
  fclose(fp);
}
  