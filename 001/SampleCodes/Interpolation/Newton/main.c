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

double  data[N+1][2], f[NUM_PNT][2], p[NUM_PNT][2];

int      winHeight=700, winWidth=700;
double   f_range, p_range;

/*------------------------------------------------
 * Test function y = cos(x)
 */
void gen_cos_data(void)
{
    int    i;

	for(i=0;i<=N;i++){
		data[i][0] = i*h;
		data[i][1] = cos(i*h*2*PI*F);
	}
}




/*--------------------------------------------------
 * Procedure to evaluate the cos() function and keep
 * NUM_PNT points in f[] array.
 */
void eval_cos(void)
{
	double  t, xmin, xmax;
    double  dx;
	int    i;

	xmin = A; xmax = B;
    dx = (xmax-xmin)/NUM_PNT;

	for(i=0;i<NUM_PNT;i++){
		t = xmin + i*dx;
		f[i][0] = t;
		f[i][1] = cos(t*2.0*PI*F);
	}
}



/*--------------------------------------------------------------------
 * Procedure to evaluate the reconstructed polynomial
 */
void eval_reconstructed(void)
{
	int      i;
	double   dx, t, xmin, xmax;

	xmin = 0.0; xmax = N*h;
	dx = (xmax-xmin)/NUM_PNT;
	t = xmin;
	for(i=0;i<NUM_PNT;i++){
		t = xmin + i*dx;
		p[i][0] = t;
		p[i][1] = Newton_interpolate(data, t, N);
	}
}


/*----------------------------------------------------------------
 * Procedure to compute the range of the real and 
 * reconstructed function values. These values are used for setting
 * up projection matrix.
 */
void compute_ranges(void)
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

	// Find the max & min reconstructed function values/
	for(i=0;i<NUM_PNT;i++){
		if(p[i][1]<pmin) pmin = p[i][1];
		if(p[i][1]>pmax) pmax = p[i][1];
	}
	p_range = pmax - pmin;
}



/*----------------------------------------------------------------------
 * The main procedure of Lagrange polynomial interpolation.
 */
int main(int argc, char **argv)
{

	// Create test data and the values of the real function.
	/*  Use cos() of different frequency, change F & N in definition.h
	    to demonstrate the sample theory and results.
    */
	gen_cos_data();    
	eval_cos();   // Compute the real cos() values.
 
	// Compute the coefficients of Newton's polynomial.
    comp_Newton_coef(data, N);
	// Evaluate the reconstructed function.
	eval_reconstructed();

/* Initialize mode and open a window in upper left corner of screen */
  glutInit(&argc,argv);             /* Make a connection to window system */
  glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB); /* Set display mode */

  /*----set window attribute ------*/
  glutInitWindowSize(winWidth,winHeight);      /* window size */
  glutInitWindowPosition(0,0);      /* initial position */
  //Initialize view parameters.
  compute_ranges();
  if(p_range<f_range)
    compute_bound_box(p, NUM_PNT);
  else
	compute_bound_box(f, NUM_PNT);

  /*----create the window ---*/
  glutCreateWindow("NewtonPolynomial");

  /*----Associate callback functions with events--*/
  glutDisplayFunc(display);         /* display event */
  glutReshapeFunc(my_reshape);      /* reshape event */
  glutKeyboardFunc(my_keyboard);    /* keyboad event */
  glutSpecialFunc(special_key_func);

  /*----Enter an infinite loop, and wait for events---*/
  glutMainLoop();
}