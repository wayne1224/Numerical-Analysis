/**************************************************************
 * Perform forward elimination for a linear system Ax=b.
 * Partial pivoting is adopted.
 */
void gauss_elm(double **A, double *b, int n);


/*----------------------------------------------------------------
 * Procedure to solve a lower triangular system
 *    L*x = b.
 */
void back_substitute(double **U, double *x, double *b, int n);
