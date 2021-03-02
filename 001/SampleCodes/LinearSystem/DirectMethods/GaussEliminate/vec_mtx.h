

/*-----------------------------------------------------
 * Procedure to allocate space for an n-dimensional
 * vector. Return the pointer to the vector.
 */
double *alloc_vec(int n);

/*---------------------------------------------------
 * Procedure to initialize an n-dimensional vector.
 */
void init_vec(double *y, int n);

/*-------------------------------------------------------------
 * Normalize a vector. If the vector norm is 0, do nothing.
 */
void normalize_vec(double *a, int n);

/*--------------------------------------------------------
 * Procedure to compute the 2-norm of a vector.
 */
double vec_norm(double *a, int n);

/*-----------------------------------------------------
 * Copy a vector from source to destination.
 *    *src: source 
 *    *dst: destination.
 */
void copy_vec(double *src, double *dst, int n);

/*--------------------------------------------------
 * Compute inner product <a, b>, where a and b are 
 * n-dimensional vectors.
 */
double  inner_product(double *a, double *b, int n);

/*----------------------------------------------------------
 * Procedure to compute the outer product of two vector.
 *   A = x*y^T, where A is an NxN matrix, and x and y are
 * n dimensional vectors.
 */
void outer_product(double **A, double *x, double *y, int n);


/*--------------------------------------------------------
 * Procedure to print out a vector.
 */
void print_vec(double *x, int n);

/*--------------------------------------------------------
 * Procedure to allocate an n X n matrix and return the 
 * pointer to the matrix.
 */
double **alloc_mtx(int n);

/*---------------------------------------------------------
 * Compute a = A*b, where a and b are vectors and A an n x n
 * matrix.
 */
void mtx_vec_mult(double *a, double **A, double *b, int n);


/*---------------------------------------------------------
 * Compute a = A^t*b, where a and b are vectors and A an n x n
 * matrix.
 */
void transpose_mtx_vec_mult(double *a, double **A, double *b, int n);

/*-------------------------------------------------------------
 * Compute B = c*A, where B and A are NxN matrices and c is 
 * a scalar.
 */
void scalar_mtx_mult(double **B, double c, double **A, int n);

/*--------------------------------------------------------
 * Procedure to compute matrix-matrix multiplication.
 *   A = B + C.
 */
void mtx_mtx_add(double **A, double **B, double **C, int n);

/*------------------------------------------------------
 * Procedure to print out a matrix.
 */
void print_mtx(double **A, int n);

/*------------------------------------------------------------
 * Procedure to make a matrix identity.
 */
void make_identity_mtx(double **A, int n);

/*-------------------------------------------------------------
 * Procedure to compute residual vector  r[] = x[] - u*y[].
 *    u: the eigen value.
 */
void comp_residual(double *r, double *x, double u, double *y, int n);



