/*************************************************************************
 * Numerical Analysis Programming #1 : 4th order Runge-Kutta Method
 * -----------------------------------------------------------------------
 * PROBLEM : Solve given differential equation by using 4th order
 *           Runger-Kutta Method with step size 0.1 and 0.09.
 *           And interpret which method is more accurate.
 *
 * INPUT   : Initial value of given differential equation
 * OUTPUT  : Numerical Solution, Exact Solution, Error 
 *
 * VARIABLE DECLARATIONS 
 *     1. double h : step size
 *     2. double t : mesh points over the interval [-1,1] 
 *     3. double x : numerical solution at t
 * -----------------------------------------------------------------------
 * Source code File: ~kkjin/Runge_Kutta/RK4/rk4.c
 * Executable  File: ~kkjin/Runge_Kutta/RK4/rk4
 * -----------------------------------------------------------------------
 * Dept. of Mathematics, Kyoungjin Kim      (Email: kkjin@math.ewha.ac.kr)
 *************************************************************************/
 #include  <stdio.h>
 #include  <math.h>

 #define   SUCCESS   1       /* given condition is true       */
 #define   LOWER_B  -1.0     /* lowerbound of interval [-1,1] */
 #define   UPPER_B   1.0     /* upperbound of interval [-1,1] */

 void      Instructions(void);             /* Function Prototype Declarations */
 double    F(double,double);
 double    Exact_Solution(double);
 void      Runge_Kutta(double,double,double);
 void      Print_Table(void);

 main()
 {
        double h,x,t;                      /* Variable Declarations */
        char   an;
        int    i;

        t=LOWER_B; x=1.0; h=0.1;           /* Initialization */

        Instructions();
        for(i=1;i<=2;i++){                 /* Apply RK4 Method with each step size */
            printf(" [%d] 4TH ORDER RUNGE-KUTTA METHOD WITH STEP SIZE %.4lf\n",i,h); 
            Runge_Kutta(t,x,h); 

	    if(i==1){                      /* pause to run RK4 with step size 0.09 */
               printf(" Press ENTER to continue..\n");
               an = getchar();
	       h-=0.01;        /* modify h to apply RK4 Method with step size 0.09 */
	    }
	    else break;
        }
 }

 /* Function : void Instructions()
  * ------------------------------
  *  Print out instructions to the user.
  */
 void   Instructions(void)
 {
        char   an;

        printf(" --------------------------------------------------------------------\n"
 	       "                    4TH ORDER RUNGE-KUTTA METHOD                     \n"
	       " --------------------------------------------------------------------\n"
	       "      We will solve given differential equation by using 4th order   \n"
	       "  Runge-Kutta with step size 0.1 and 0.09 over the interval [-1,1]   \n"
	       "                                                                     \n"
	       "                        x' = x+t, if -1<=t<=0                        \n" 
	       "                             x-t, if  0<=t<=1                        \n"
	       "                        x(-1) = 1 (initial value)                    \n");
        printf(" --------------------------------------------------------------------\n");
        printf(" Press ENTER to continue..\n"); 
        an = getchar();
 }


/* Function : double F(double,double)
 * -----------------------------------------
 *  Calculate function value at (t,x) and return the result. 
 */
 double F(double t,double x)
 {
        if(t==0.0)      return  x;
        else if(t<0.0)  return (x+t);
        else            return (x-t);
 }

/* Function : double Exact_Solution(double)
 * ----------------------------------------
 *  Calculate exact solution of given differential equation at t.
 */
 double Exact_Solution(double t)
 {
        if(t<=0.0) return (exp(t+1)-t-1);
        else       return (exp(t+1)-2.0*exp(t)+t+1);
 }

/* Function : void Runge_Kutta(double,double,double)
 * -------------------------------------------------
 *  Find a solution of given differential equation using 4th order Runge-Kutta Method. 
 *         
 *  RUNGE-KUTTA METHOD OF ORDER 4
 *  -----------------------------              
 *              x(t+h) = x(t)+(F1+2F3+2F3+F4)/6 
 *
 *                  F1 = hf(t,x)
 *                  F2 = hf(t+h/2, x+F1/2)
 *                  F3 = hf(t+h/2, x+F2/2)
 *                  F4 = hf(t+h,x+F3) 
 */
 void   Runge_Kutta(double t,double x,double h)
 {
        static double f[4];         /* array to store F1 to F4 */
        double error,exact,e_sum;
        int    i=1;

        exact=Exact_Solution(t);    /* At zero step : initial value */
        error=fabs(x-exact);
        e_sum=error;
        Print_Table(); 
        printf("%4d %15.8lf %15.8lf %15.8lf %15.8lf\n",0,t,x,exact,error);

        while(SUCCESS){             /* 4th order Runge-Kutta Method */
           f[0] = h*F(t,x);
           f[1] = h*F(t+h/2.0,x+f[0]/2.0);
           f[2] = h*F(t+h/2.0,x+f[1]/2.0);
           f[3] = h*F(t+h,x+f[2]);
	   
	   x += (f[0]+2*(f[1]+f[2])+f[3])/6.0;
           t = LOWER_B + h*i;                 /* update t */

	   exact = Exact_Solution(t);             /* find exact solution at t */
	   error = fabs(x-exact);                 /* calculate error */
           e_sum+=error;                          /* update error sum*/
           
	   if(t<=UPPER_B)             /* print solution and update step number*/
              printf("%4d %15.8lf %15.8lf %15.8f %15.8lf\n",i++,t,x,exact,error);
	   else break;                /* t exceeds given range -> exit */
        }
        printf(" --------------------------------------------------------------------\n");
        printf(" Total truncation & roundoff error : %13.8lf\n",e_sum);
        printf(" Average of local error            : %13.8lf\n\n",e_sum/i);
 }

/* Function : void Print_Table()
 * -----------------------------
 *  Print the header of the solution table.
 */
 void   Print_Table(void)
 {
        printf(" ====================================================================\n");
        printf("   n          t              x(t)          Real Sol          Error\n");
        printf(" --------------------------------------------------------------------\n");
 }
