/****************************************************************************
 * This sample program reveals the accuracy of Gaussian quadrature integration method.
 * The method: different numbers of intervals X different numbers of Gaussian points
 */
#include <stdio.h>
#include <math.h>

//max number of interval
#define   N  10

//max number of Gaussian points
#define  K   4

//Sample points in [-1, 1], 1-, 2-, 3-, and 4-points. Ignore extra points (the 0s).
double  p[K][K]={{0.0, 0.0, 0.0, 0.0}, {0.5773502691896257, -0.5773502691896257, 0.0, 0.0}, 
                           {0.0, 0.7745966692414834, -0.7745966692414834, 0.0}, 
                           {0.3399810435848563, -0.3399810435848563, 0.8611363115940526, -0.8611363115940526}};
//function values in N interval
double  f[N][K];
//weights, for 1-, 2-, 3-, and 4-points. Extra weights = 0s.
double  wgt[K][K]={{2.0, 0.0, 0.0, 0.0}, {1.0, 1.0, 0.0, 0.0}, 
                                 {0.888888888888888888888, 0.555555555555555556, 0.555555555555555556, 0.0},
                                 {0.6521451548625461, 0.6521451548625461, 0.3478548451374538, 0.3478548451374538}};
//integral range = [0, 1].
double  a=0.0, b=3.0;


/*------------------------------------------------------------------------------------------------------
 * The function to be integrated: f(x) = 1.0/(1+x^2).
 * The indefinite integral I(f(x)=F(x) =  arctan(x)
 */
double  my_func(double x)
{
   double  y;
   int         k;

   /*
   y = 1.0;
   for(k=4;k>=0;k--)
      y = y*x + 1.0;
   return(y);
   */
  return(2.0*x/(1.0+x*x));
}


/*-----------------------------------------------------------------------------------------------------
 * The inetgration of f(x) = 1.0/(1+x^2).
 *  F(x) =arctan(x)
 */
double my_func_integral(double  a, double b)
{
   double     t;
   double     Fa, Fb;
   int            i;

   Fa = a/6.0;
   Fb = b/6.0;
   for(i=5;i>=1;i--){
	   Fa = Fa*a + a/i;
	   Fb = Fb*b + b/i;
   }
   /*
    Fa = atan(a);
     Fb =atan(b);
	 */
   return(Fb-Fa);
}
  
/*-----------------------------------------------------------------------------------------------------
 * Conformation mapping from [-1, 1] to [a, b]
 *  Input:
 *     t: Gaussian point in [-1, 1].
 *     a, b: the 2 ends of the interval.
 *  Return: 
 *     Gaussian point in [a, b].
 *  Called before evaluating function values.
 */
double  conform_map(double t, double a, double b)
{
   double   x;

   x = (b-a)*t/2.0 + (a+b)/2.0;
   return (x);
}


/*---------------------------------------------------------------------------------------------------------
 * Gaussian quadrature integration subroutine.
 * Input:
 *   k: number of Gaussian points used,
 *   a, b: 2 ends of the interval,
 *   f(): the funtion.
 * Return:
 *   The integral value.
 */
double  gauss_quadrature(int k, double a, double b, double f())
{
    int   i, j;
    double   sum=0.0;
    double   x, t;

    //Acumulate the contribution of all Gaussian points.
    for(i=0;i<k;i++){
        t = p[k-1][i]; //Retrieve the cononical Gauss point
        x = conform_map(t, a, b); //Conformation mapping
        sum += wgt[k-1][i]*f(x); //The integration
    }
     sum = sum*(b-a)/2.0; //Scale the integration value.
    return (sum);
}

/*----------------------------------------------------------------------------------------------------------
 * The main procedure 
 */
void main(int argc, char **argv)
{
  char      filename[]="gauss.txt";
  int        k, numInterval, i, j;
  double  gaussSum, localSum;
  double   exactSum, h;
  FILE     *fp;


  //Open the output file
  fp = fopen(filename, "w");
  //Compute the exact sum.
 // exactSum = my_func_integral(a, b);
//  fprintf(fp,"Exact integration = %15.10f\n", exactSum);
  fprintf(fp,"-----------------------------------------------------------------------------\n");
  //Using different number of Gaussian points
  for(k=1;k<=4;k++){
	  fprintf(fp, "No. of Gaussian points = %d\n", k);
  //Using different numbers of intervals
      for(numInterval=1; numInterval<=N; numInterval++){
	      fprintf(fp, "  No. of interval = %d\n", numInterval);
          h = (b-a)/numInterval;
	     gaussSum = 0.0;
		  //Integrate function in each interval
	      for(i=0;i<numInterval;i++){
              localSum = gauss_quadrature(k, a+i*h, a+(i+1)*h, my_func);
		      gaussSum += localSum;
		  }
          fprintf(fp,"    Integral= %15.10f\n",gaussSum);
	  }//end for(numInterval)
  }//end for(k)
  fclose(fp);
}//end of main()
  