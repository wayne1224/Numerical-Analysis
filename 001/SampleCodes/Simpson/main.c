/***************************************************************
 * Main procedure of the Lagrange polynomial interpolation method.
 * This program produces test data at discrete posistions and 
 * reconstruct the function. Then a graphics routine shows the 
 * real function and the reconstructed function.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <GL/glut.h>


#include "definition.h"

// Sample points.
double  data[MAX][2];
// Dense sample points for drawing function.
double  f[NUM_PNT][2], f_range, p[NUM_PNT][2];

int      winHeight=700, winWidth=700;
int      N;   // last index of sample point.

double   exactVal,  // exact solution.
         simpsonVal;  // solution of simpson's rule.
double   h;           // sample interval size.
char     mystring[64];

/*------------------------------------------------
 * Test function y = cos(f*x) + 1.0
 *   frequency: F.
 *   A<= x <=B.
 */
void gen_sample_data(void)
{
    int    i;

	for(i=0;i<=N;i++){
		data[i][0] = i*h;
		data[i][1] = cos(F*i*h) + 1.0;
	}
}




/*--------------------------------------------------
 * Procedure to evaluate the cos(F*x) function and keep
 * NUM_PNT points in f[] array.
 */
void eval_func_values(void)
{
	double  t, xmin, xmax;
    double  dx;
	int     i;

	xmin = A; xmax = B;
    dx = (xmax-xmin)/NUM_PNT;
    t = xmin;
	for(i=0;i<NUM_PNT;i++){
		f[i][0] = t;
		f[i][1] = cos(F*t) + 1.0;
		t += dx;
	}
}


/*-----------------------------------------------------
 * Procedure to evaluate the interpolation polynomials of
 * Simpson's rule.
 *   data[][2]: the sampled data,
 *   n: the last index of sample point.
 * NUM_PNT points in p[] array.
 */
void eval_polynomial(double data[][2], int n)
{
	double  t, xmin, xmax;
	double  a, b, c, fa, fb, fc, func;
    double  dx;
	int     i, j;

	xmin = A; xmax = B;
    dx = (xmax-xmin)/NUM_PNT;
    t = xmin;
	j = 0;
	for(i=0;i<n;i+=2){
		a = data[i][0]; b= data[i+1][0]; c=data[i+2][0];
		fa = data[i][1]; fb = data[i+1][1];
		fc = data[i+2][1];
        while(t<c){
           func = fa*(t-b)*(t-c)/(2.0*h*h) +
			      fb*(t-a)*(t-c)/(-h*h) +
				  fc*(t-a)*(t-b)/(2.0*h*h);
		   p[j][0] = t;
		   p[j][1] = func;
		   j++;
		   t = t + dx;
		}
	}
}

/*--------------------------------------------------------------------
 * Procedure to evaluate the exact solution of the integration.
 *    g(x) = cos(2fx)+1.
 *    Int(g(x)) = sin(fx)/(f) + x.
 */
double eval_exact_integral(void)
{
	double   GA, GB;

    GA = sin(F*A)/(F) + A;
    GB = sin(F*B)/(F) + B;
	return(GB-GA);
}


/*----------------------------------------------------------------
 * Procedure to compute the range of the  sample points in 
 * y directions. The ranges are used to setting up the projection
 * matrix.
 */
void compute_func_range(void)
{
	double  fmin=1.0e30, fmax=-1.0e30;
    double  pmin=1.0e30, pmax=-1.0e30;
	int     i;

    // Find the max and min function values.
	for(i=0;i<NUM_PNT;i++){
		if(f[i][1]<fmin) fmin = f[i][1];
		if(f[i][1]>fmax) fmax = f[i][1];
	}
	f_range = fmax - fmin;
}



/*----------------------------------------------------------------------
 * The main procedure of Lagrange polynomial interpolation.
 */
int main(int argc, char **argv)
{
    
	// Create test data.
	// Compute the values of the function.
	/*  Use cos() of different frequency, change F in definition.h
	    to demonstrate the sample theory and results.
    */
	N = 4;
	h = (B-A)/N;
	gen_sample_data();    
	eval_func_values();   // Compute the real function values.
    compute_func_range();

	// Evaluate the interpolation function.
    eval_polynomial(data, N);

	// Evaluate the reconstructed function.
	exactVal = eval_exact_integral();

	// Carry out trapezoid rule to integrate the function.

	simpsonVal = simpson(data, h, N);
    sprintf(mystring,"N=%d, h=%lf, exact Val=%lf, simpson=%lf, err=%lf\0", 
		N, h, exactVal, simpsonVal, fabs(exactVal-simpsonVal));
/* Initialize mode and open a window in upper left corner of screen */
  glutInit(&argc,argv);             /* Make a connection to window system */
  glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB); /* Set display mode */

  /*----set window attribute ------*/
  glutInitWindowSize(winWidth,winHeight);      /* window size */
  glutInitWindowPosition(0,0);      /* initial position */

  //Initialize view parameters.
  compute_bound_box(f, NUM_PNT);

  /*----create the window ---*/
  glutCreateWindow("Simpson's Rule");

  /*----Associate callback functions with events--*/
  glutDisplayFunc(display);         /* display event */
  glutReshapeFunc(my_reshape);      /* reshape event */
  glutKeyboardFunc(my_keyboard);    /* keyboad event */
  glutSpecialFunc(special_key_func);

  /*----Enter an infinite loop, and wait for events---*/
  glutMainLoop();
}