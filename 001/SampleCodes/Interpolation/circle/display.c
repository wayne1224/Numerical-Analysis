/****************************************************************
 * Procedure to draw a list of 2D points in the xy-planes.
  */
#include <stdio.h>

#include <windows.h>
#include <GL/glut.h>

#include "mydefinition.h"

float   left, right, top, bottom;
float   xcenter, ycenter;
float   eye[3], focus[3];
float   xmin, ymin, xmax, ymax;

// External variables declared in main.c

extern int  winHeight, winWidth;
extern double  data[N+1][2], f[NUM_PNT][2], p[NUM_PNT][2];

/*----------------------------------------------------------------------
 * Procedure to compute the bounding box of the function graph
 */
void compute_bound_box(double  pnt[][2], int size)
{  
   int   i;



   ymin = xmin = 1.0e30;
   ymax = xmax = -1.0e30;
   for(i=0;i<size;i++){
      if(xmin > pnt[i][0]) xmin = pnt[i][0];
      if(ymin > pnt[i][1]) ymin = pnt[i][1];
      if(xmax<pnt[i][0]) xmax = pnt[i][0];
      if(ymax<pnt[i][1]) ymax = pnt[i][1];
   }
   
   fprintf(stderr,"bounding box= %lf %lf --> %lf %lf\n", xmin, ymin, xmax, ymax);
   focus[0] = eye[0] = xcenter = (xmin + xmax)/2.0;
   focus[1] = eye[1] = ycenter = (ymin + ymax)/2.0;
   left = xmin - xcenter; right = xmax - xcenter;
   bottom = ymin - ycenter; top = ymax - ycenter;
}


/*------------------------------------------------------
 * Procedure to draw points.
 */
void draw_func(double pnt[][2], int size)
{
	int  i;


	for(i=0;i<size-1;i++){
      glBegin(GL_LINES);
	    glVertex3f(pnt[i][0], pnt[i][1], 0.0);
		glVertex3f(pnt[i+1][0], pnt[i+1][1], 0.0);
	  glEnd();
	}

}


/*-----------------------------------------------------------
 * Procedure to draw sample points.
 *   n: last index of the sample points.
 */
void draw_sample_points(double pnt[][2], int n)
{
	int  i;


	glColor3f(0, 1, 0);
	glPointSize(5.3);
	glBegin(GL_POINTS);
	  for(i=0;i<=n;i++) glVertex3f(pnt[i][0], pnt[i][1], 0.0);
	glEnd();
}


/*---------------------------------------------------
 * Procedure to draw grid.
 */
void draw_grid_lines(void)
{
   int i;
   float dx, dy;
   float x, y;

   dx = (xmax-xmin)/NUM_LINE;
   dy = (ymax-ymin)/NUM_LINE;

   glColor3f(0.3, 0.3, 0.3);
   x = xmin; y = ymin;
   for(i=0;i<=NUM_LINE;i++){
	   glBegin(GL_LINES);
	     glVertex3f(x, ymin, 0);
		 glVertex3f(x, ymax, 0);
	   glEnd();
	   x += dx;

	   glBegin(GL_LINES);
	     glVertex3f(xmin, y, 0);
		 glVertex3f(xmax, y, 0);
	   glEnd();
	   y += dy;
   }
   //Draw the x & y axes.
   glColor3f(0, 0, 0);
   glLineWidth(3.4);
   glBegin(GL_LINES);
	 glVertex3f(-100, 0, 0);
	 glVertex3f(100, 0, 0);
   glEnd();
   glBegin(GL_LINES);
	 glVertex3f(0, -100, 0);
	 glVertex3f(0, 100, 0);
   glEnd();

}

/*------------------------------------------------------
 * Procedure to draw the output string.
 */
void draw_string(char *mystring, int x, int y, float r, float g, float b)
{
	int i, len;

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
    glOrtho(0.0, winWidth, 0.0, winHeight, -1.0, 1.0);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	len = strlen(mystring);
	glColor3f(r, g, b);
	for(i=0;i<len;i++){
        glRasterPos2i(x+9*i, y);
		glutBitmapCharacter(GLUT_BITMAP_8_BY_13, mystring[i]);
	}
	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
}

/*-------------------------------------------------
 * Callback function of reshape event
 *   w: new width of the window, 
 *   h: new height of the window
 */
void my_reshape(int width, int height)
{
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
 /*---- Record the dimension of the window ---*/ 
  winWidth = width;
  winHeight = height;

  glOrtho(left, right, bottom, top, -1.0, 1.0);
 // glOrtho(0, xmax, 0, ymax, -10.0, 10.0);

  glViewport(0, 0, width, height);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glutPostRedisplay();     /* trigger a display event */
}


/*----------------------------------------------
 * Callback function for keyboard event.
 */
void my_keyboard(unsigned char key, int x, int y)
{

  if(key=='q') exit(0);
  else if(key=='n'){
    left = 0.95*left;
    right = 0.95*right;
    top = 0.95*top;
    bottom = 0.95*bottom;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(left, right, bottom, top, -100, 100); 
  }else if(key=='f'){
    left = 1.05*left;
    right = 1.05*right;
    top = 1.05*top;
    bottom = 1.05*bottom;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(left, right, bottom, top, -100, 100); 
  }
  display();
}

/*--------------------------------------------------
 * Callback for special key
 */
void special_key_func(int key, int x, int y)
{
  float  dx, dy;

  dx = (right-left)/20.0;
  dy = (top-bottom)/20.0;
  if(key==GLUT_KEY_LEFT){
    eye[0] += 0.3*dx;
    focus[0] += 0.3*dx;
  }else if(key==GLUT_KEY_RIGHT){
    eye[0] -= 0.3*dx;
    focus[0] -= 0.3*dx;
  }else if(key==GLUT_KEY_UP){
    eye[1] += 0.3*dy;
    focus[1] += 0.3*dy;
  }else if(key==GLUT_KEY_DOWN){
    eye[1] -= 0.3*dy;
    focus[1] -= 0.3*dy;
  }
  display();
}


/*-----------------------------------------------------------------------------------------
 * Callback of the display event.
 * It clears the window and draws the sample points and the interpolation results.
 */
void my_display(void)
{ 
    glClearColor(0.9, 0.9, 0.9, 0.0);  //white background
    glClear(GL_COLOR_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(eye[0], eye[1], 1.0, focus[0], focus[1], 0.0, 0, 1, 0);

	draw_grid_lines();
	glColor3f(1,0, 0);
	draw_func(f, NUM_PNT);
	glColor3f(0, 0, 1);
	draw_func(p, NUM_PNT);
    draw_sample_points(data, N);

	draw_string("Actual func.", winWidth/2-10, 20, 1, 0, 0);
	draw_string("Interpolation func.", winWidth/2-10, 35, 0, 0, 1);
	glutSwapBuffers();
}