//Define PI
#define  PI 3.141592653

// Define the 2 ends of the domain.
#define  A  0.0
#define  B  1.6*PI

#define  MAX  256

// Define the number of points generated when drawing the functions.
#define  NUM_PNT    500

//define frequency of cos(x)
#define  F  3

/*------------------------------------------------
 * Test function y = cos(f*x) + 1.0
 *   frequency: F.
 *   A<= x <=B.
 */
void gen_sample_data(void);

/*-----------------------------------------------------
 * Procedure to evaluate the interpolation polynomials of
 * Simpson's rule.
 *   data[][2]: the sampled data,
 *   n: the last index of sample point.
 * NUM_PNT points in p[] array.
 */
void eval_polynomial(double data[][2], int n);


/*---------------------------------------------------------
 * Procedure to compute the integral by using trapezoid rule.
 *   data[][2]: sample points, (xi, yi)
 *   h: gap between sample points,
 *   n: largest index of sample point.
 * Return the integral value.
 */
double simpson(double data[][2], double h, int n);


/*----------------------------------------------------------------------
 * Procedure to compute the bounding box of the function graph
 */
void compute_bound_box(double  pnt[][2], int size);

/*-------------------------------------------------
 * Callback function of reshape event
 *   w: new width of the window, 
 *   h: new height of the window
 */
void my_reshape(int width, int height);

/*----------------------------------------------
 * Callback function for keyboard event.
 */
void my_keyboard(unsigned char key, int x, int y);

/*--------------------------------------------------
 * Callback for special key
 */
void special_key_func(int key, int x, int y);


/*--------------------------------------------------
* display event callback function.
*/
void display(void);

