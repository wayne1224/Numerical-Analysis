//Define PI
#define  PI 3.14159

// Define the 2 ends of the domain.
#define  A  0.0
#define  B  1.0

// define the max index of sample points.
// There are (N+1) sample points.
#define  N  52

// Define the number of points generated when drawing the functions.
#define  NUM_PNT    400

//define gap between sample points.
#define  h          (B-A)/N




//define frequency of cos(x)
#define  F 4

/*---------------------------------------------------------
 * Procedure to compute the coefficients of Newton
 * polynomial.
 *   data[][2]: sample points, (xi, yi)
 *   n: largest index of sample point.
 * Return the interpolated value.
 */
void comp_Newton_coef(double data[][2], int n);

/*-------------------------------------------------------------
 * Procedure to perform Newton polynomial interpolation.
 *   data[][2]: the sample points,
 *   t: the input parameter value,
 *   n: last index of sample pnt.
 * Return the interpolayed value.
 */
double Newton_interpolate(double data[][2], double t, int n);


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

