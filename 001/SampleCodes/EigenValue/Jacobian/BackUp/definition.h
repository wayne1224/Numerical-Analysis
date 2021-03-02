#define   EPSILON     0.000001
#define   MAXSTEP     100

#define   PI          3.141592653


/*-------------------------------------------------------------
 * Jacobian method to compute the eigenvalues and eigenvectors.
 *   A: the input matrix,an n X n mtx,
 *   P: matrix of the eigen vectors,
 *   n: dimension of the matrix.
 * Matrix A will be diagonalized. The eigen values are at the
 * main diagonal.
 */
void jacobian_method(double  **A, double **P,  int n);


