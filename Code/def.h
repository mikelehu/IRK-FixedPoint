/*----------------------------------------------------------------------------*/
/*									      */
/*                                def.h                                       */
/*					                                      */
/* ---------------------------------------------------------------------------*/

#include "prec.h"
#include <stdbool.h>
#include <float.h> 

/* ---------------------------------------------------------------------------*/
/*									      */
/*	Parameters							      */
/*									      */
/* ---------------------------------------------------------------------------*/

#ifndef DEFH_
#define DEFH_

#define PI atan(1.) * 4.   
#define INF DBL_MAX
#define MAXNEQ	200		// Maximum neq.
#define MAXPARAM 30		// Number of maximum parameters.
#define MAXIT 100		// Maximum number of fixed point iterations 
#define MAXKSW 10               // Max number of steps to change second 
                                // integration initialization mode.
#define SUCCESS 0
#define FAIL -1
#define STRMAX 40		// Filename maximum string.
#define RTOL pow(10.,-12);	// pow(2.,-40);
#define ATOL pow(10.,-12);
#define ESTERRTHRS pow(10.,16)	// EstimatedErrorThreshold. pow(10.,12)
#define IOUT
#define PARALLEL		// Active opem-mpi parallel execution 

#define DIR_TERM "../CoefficientsData/"    // Path Coefficients (terminal)
#define DIR_MATH "../../CoefficientsData/" // Path Coefficients (mathematica)

#define DOUBLEPRECISION
typedef double val_type;


//DOUBLEPRECISION or FLOAT

#define POW(x, y)      pow(x, y)
#define SQRT(x)        sqrt(x)

#define EXP(x)        exp(x)
#define LOG(x)        log(x)

#define SIN(x)        sin(x)
#define COS(x)        cos(x)
#define TAN(x)        tan(x)

#define FABS(x)	      fabs(x)
#define FMAX(x,y)     fmax(x,y)

/* ---------------------------------------------------------------------------*/
/*									      */
/*	General definitions					       	      */
/*								              */
/* ---------------------------------------------------------------------------*/

typedef struct gauss_method
   {
     int ns;
     val_type *m;	 	 // mij=aij/bj and mij+mji-1=0.
     val_type *c,*b,*a;	 	 // c,b,a coefficients.
     val_type *hc;       	 // hc=h*c.
     val_type *hb;       	 // hb=h*b.
     val_type *nu;               // interpolate coeficcients (a*/bj).
     int *orderedindices;	 // ascending order indices for bi coeficients

   } gauss_method;


typedef struct solution
  { 
     val_type *uu,*ee;

  } solution;


typedef struct toptions
  {
     val_type h; 	        // stepsize
     val_type t0;
     val_type t1;
     val_type *rtol,*atol;
     int algorithm;			  
     int sampling;
     int approximation;
     int rdigits,mrdigits;
     int (*iteration[2])();     // Iteration : Jacobi or Seidel.
   } toptions; 


typedef struct parameters
   { 
     int numrpar;	 // Number of real parameters.
     val_type *rpar;	 // Variables for specifying odefun real parameters. 
     int numipar;	 // Number of int parametes. 
     int *ipar;          // Variables for specifying odefun integer parameters.
     int eval;	         // specify which part of differential equation 
                         // must be evaluated

   } parameters;


typedef struct ode_sys       
   {
     int problem;
     int neq;		// number of equations.
     int n;		// n-body.
     void (*f)();	// odefun.
     val_type (*ham)();	// hamiltonian
     int cod[2];	// tO specify which part of ODE must be evaluated.
     parameters params;

    } ode_sys;


typedef struct solver_stat
    {

    int convergence;	      // SUCCESS or FAIL. 			        
    bool laststep;                                

    /* auxiliar*/
    val_type *z,*li,*zit0,*fz;

    val_type E0;     	      // Initial Energy
    val_type MaxDE; 	      // MaxDE=Abs(Ei-E0/E0)    
		 
    /* stadistics*/
    int stepcount;
    int itcount;
    long int totitcount;
    long int totitcountzero;	      // 28-11-2016.
    int maxitcount;			
    int itzero;				
    int fcn;  
    int *initqlty;	      // ns*neq matrix. Quality of initialization of Yi stages 

    /* output filename */
    char filename[STRMAX];  // Integration filename.
    int nout;               // number of output values.

    } solver_stat; 
	
#endif /*DEFH_*/
