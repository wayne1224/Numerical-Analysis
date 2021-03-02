/*------------------------------------------------------
 * Procedure to do similarity transformation:
 *   A = PAP, 
 *   Q = Q*P.
 * Such that A becomes a Hessenberg matrix.
 */
void Hessenberg_form(double **A, double **Q, int n);

/*---------------------------------------------------------------
 * The procedure to compute the eigen-values.
 * Input:
 *    T: a tri-diagonal matrix,
 *    n: dimension of the mtx.
 * Output:
 *    eigenVal[]: array of the eigen-values,
 *    numEigen: num. of eigen-values.
 */
void tri_diag_mtx_eigen(double **T, int n, double *eigenVal, int *numEigen);