#define  EPSILON      0.000001

/*---------------------------------------------------
 * Procedure to solve a linear system by using Conjugate
 * Gradient (CG) method.
 * Algm:
 *    1. Select x0;
 *    2. Compute d = b - Ax0;
 *    3. Compute g = -d;
 *    4. For k=0, 1, 2, 3, ..., n-1
 *           h = A*d;
 *          Compute alpha;
 *          x = x + alpha*d;
 *          g = g + alpha*h;
 *          Compute beta;
 *          d = -g + beta*d;
 *   5. end-for
 */
void conjugate_gradient(double **A, double *x, double *b, int n);

