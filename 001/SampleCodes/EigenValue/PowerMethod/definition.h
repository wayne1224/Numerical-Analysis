#define   EPSILON     0.0001
#define   MAXSTEP     40

/*-------------------------------------------------------------
 * Power method to compute the max eigenvalue and eigenvector.
 *   A: an n X n mtx,
 *   x: the max eigen vector,
 *   val: the pointer to the max eigen value,
 *   n: dimension of the matrix.
 */
void power_method(double  **A, double *x, double *val, int n);

/*-------------------------------------------------------------
 * Inverse Power method to compute the min eigenvalue and eigenvector.
 *   A: an n X n mtx,
 *   x: the min eigen vector,
 *   val: the pointer to the min eigen value,
 *   n: dimension of the matrix.
 */
void inverse_power_method(double  **A, double *x, double *val, int n);

/*-------------------------------------------------------------------
 * Procedure to perform Doolittle LU decomposition.
 */
void doolittle(double **A, double **L, double **U, int n);
