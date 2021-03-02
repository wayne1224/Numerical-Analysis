/*********************************************************************
 * Kernel procedure of the conjugate gradient method.
 * Called by main() to solve Ax = b in main.c.
 */
#include <stdio.h>

#define EPSILON      0.00000001
#include "conjugate_gradient.h"
#include "vec_mtx.h"

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
void conjugate_gradient(double **A, double *x, double *b, int n)
{
   double       *h, *d, *g, *t;
   double        newG2, oldG2, dh;
   double        error, alpha, beta;
   int               k;

   //Allocate memory spaes for the temperary vectors
   h = alloc_vec(n);
   d = alloc_vec(n);
   g = alloc_vec(n);
   t = alloc_vec(n);

  //Initialize x0
   init_vec(x, n);
   
  //Compute d0 = b - Ax0;
   mtx_vec_mult(t, A, x, n);  // t = Ax.
   vec_vec_sub(d, b, t, n);     // d = b-t = b - Ax.
   scalar_vec_mult(g, -1.0, d, n);  //g = -d.
   error = vec_norm(d, n);

   fprintf(stderr,"Initial x[]=\n");
   print_vec(x, n);
   fprintf(stderr,"d= ");
   print_vec(d, n);
   fprintf(stderr,"g= ");
   print_vec(g, n);
   //The iterations
   for(k=0;error>EPSILON;k++){
      mtx_vec_mult(h, A, d, n);         // h = A*d.
      dh = inner_product(d, h, n);      // <d,h>.
      oldG2 = inner_product(g, g, n);
      alpha = oldG2/dh;

      // Update x[]
      scalar_vec_mult(t, alpha, d, n);     //t = alpha*d
      vec_vec_add(x, x, t, n);                //x = x + alpha*d
	  error = vec_norm(t, n);
      // Update gradient g[]
      scalar_vec_mult(t, alpha, h, n);     //t = alpha*h
      vec_vec_add(g, g, t, n);                //g =g + alpha*h

    // Compute new searching direction d[].
      newG2 = inner_product(g, g, n);   //newG2 = <g, g>
      beta = newG2/oldG2;
      scalar_vec_mult(t, beta, d, n);      //t = beta*d;
      vec_vec_sub(d, t, g, n);               //d = -g + beta*d

       //Print out the vector
	   fprintf(stderr,"---------------------------------------------\n");
       fprintf(stderr,"k=%d, x[]=\n", k);
       print_vec(x, n);
	   fprintf(stderr,"alpha = %lf, beta = %lf\n", alpha, beta);
       fprintf(stderr,"d= ");
       print_vec(d, n);
	   fprintf(stderr,"g= ");
       print_vec(g, n);
   }//next for(k)
}//end_proc


      
      