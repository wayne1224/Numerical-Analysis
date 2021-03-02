/*******************************************************************************
 * This file contains the main procedure of this example program. The main procedure performs
 * the following computation:
 * 1. generating the sample points,
 * 2. computing the x- and y-coordinates of NUM_INTER_PNTS interpolation points,
 * 3. displaying the results.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <GL/glut.h>

#include "mydefinition.h"

//Number of samples in a circle
int    numSamples = 12;
//Radius of the circle  
double  radius = RADIUS;

//Sample points' space. A large chunk of memory
double  x[NUM_INTER_PNTS], y[NUM_INTER_PNTS], t[NUM_INTER_PNTS];

//Interpolation points' memory
double  px[NUM_INTER_PNTS], py[NUM_INTER_PNTS];
 
/*---------------------------------------------------------------------------------------------------------
 * This is the main procedure. It works as follows:
 *   1. alocating memory space
void main(int argc, char **argv)
{
//generate the sample points.
   gen_circle_pnts(t, x, y, numSamples);
   //Compute x- and y-coordinates of NUM_INTER_PNTS interpolation points.
   Lagrange(numSamples, t, x, px);
   Lagrange(numSamples, t, y, py);
  
   //Initialize glut
   glutInit(&argc, argv);
   glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB);
   //Create the window
   glutInitWindowSize(800, 800);
   glutInitWindowPosition(0, 0);
   glutCreateWindow("Circle using Lagrange");
   //Register callbacks
   glutDisplayFunc(my_display);
   glutReshapeFunc(my_reshape);
   glutKeyboardFunc(my_keyboard);

    //Enter the arbiter.
    glutMainLoop();

}