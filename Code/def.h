/*------------------------------------------------------------------------------*/
/*										*/
/*                                def.h    					*/
/*										*/
/* -----------------------------------------------------------------------------*/
#include "prec.h"
#include <stdbool.h>

/* -----------------------------------------------------------------------------*/
/*										*/
/*	Parameters								*/
/*										*/
/* -----------------------------------------------------------------------------*/

#ifndef DEFH_
#define DEFH_

#define PI atan(1.) * 4.   
#define INF 99999.
#define MAXNEQ	200				// Maximum neq.
#define MAXPARAM 30				// Number of maximum parameters.
#define MAXIT 50				// Maximum number of fixed point iterations 50
#define MAXKSW 10                               // Max number of steps to change second 
                                                // integration initialization mode.
#define SUCCESS 0
#define FAIL -1
#define STRMAX 40				// Filename maximum string.
#define RTOL pow(10.,-12);			// pow(2.,-40);
#define ATOL pow(10.,-12);
#define ESTERRTHRS pow(10.,16)			// EstimatedErrorThreshold. pow(10.,12)
#define IOUT
#define PARALLEL				// Active opem-mpi parallel execution     

/* -----------------------------------------------------------------------------*/
/*										*/
/*	PREC: FLOAT,DOUBLEPRECISION,QUADRUPLEPRECISION				*/
/*										*/
/* -----------------------------------------------------------------------------*/

//#define PREC 2        to be specified in the Makefile file ( -D PREC=value)
//1=DOUBLEPRECISION    				// Specify by: gcc -D DOUBLE ...
//2=QUADRUPLEPRECISION			  
//3=FLOAT

#if PREC ==2 
#define QUADRUPLEPRECISION 
typedef __float128 val_type;
typedef double lowfloat;
#elif PREC ==3
#define FLOAT
typedef float val_type;
typedef float lowfloat;
#define RTOL pow(2.,-16);
#define ATOL pow(2.,-16);
#else
#define DOUBLEPRECISION
typedef double val_type;
typedef float lowfloat;
#endif


#ifdef QUADRUPLEPRECISION

#define POW(x, y)     powq(x, y)
#define SQRT(x)       sqrtq(x)

#define EXP(x)        expq(x)
#define LOG(x)        logq(x)

#define SIN(x)        sinq(x)
#define COS(x)        cosq(x)
#define TAN(x)        tanq(x)

#define FABS(x)	      fabsq(x)
#define FMAX(x,y)     fmaxq(x,y)

#else //DOUBLEPRECISION or FLOAT

#define POW(x, y)      pow(x, y)
#define SQRT(x)        sqrt(x)

#define EXP(x)        exp(x)
#define LOG(x)        log(x)

#define SIN(x)        sin(x)
#define COS(x)        cos(x)
#define TAN(x)        tan(x)

#define FABS(x)	      fabs(x)
#define FMAX(x,y)     fmax(x,y)

#endif

/* -----------------------------------------------------------------------------*/
/*										*/
/*	General definitions							*/
/*										*/
/* -----------------------------------------------------------------------------*/

typedef struct gauss_method
   {
     int ns;
     val_type *m;	 			 // mij=aij/bj and mij+mji-1=0.
     val_type *c,*b;	 			 // c,b coefficients.
     val_type *hc;       			 // hc=h*c.
     val_type *hb;       			 // hb=h*b.
     val_type *nu;                               // interpolate coeficcients (a*/bj).
     int *orderedindices;			 // ascending order indices for bi coeficients
   } gauss_method;


typedef struct solution
  { 
     val_type *uu,*ee;

  } solution;


typedef struct toptions
  {
     val_type h; 				  // stepsize
     val_type t0;
     val_type t1;
     val_type *rtol,*atol;
     int algorithm;			  
     int sampling;
     int approximation;
     int rdigits,mrdigits;
     int (*iteration[2])();			 // Iteration : Jacobi1, Jacobi2.
   } toptions; 


typedef struct parameters
   { 
     val_type *rpar;				 // Variables for specifying odefun real parameters. 
     int *ipar;     			         // Variables for specifying odefun integer parameters.
						 //  ipar[0]: which part of differential equation must be evaluated

   } parameters;


typedef struct ode_sys       
   {
     int problem;
     int neq;					// number of equations.
     int n;					// n-body.
     void (*f)();				// odefun.
     val_type (*ham)();				// hamiltonian
     int cod[2];				// tO specify which part of ODE must be evaluated.
     parameters params;

    } ode_sys;


typedef struct solver_stat
    {

    int convergence;				  // SUCCESS or FAIL. 			        
    bool laststep;                                

    /* auxiliar*/
    val_type *z,*li;

    val_type E0;     		      		  // Initial Energy
    val_type MaxDE;		       		  // MaxDE=Abs(Ei-E0/E0)    
		 

    /* stadistics*/
    int stepcount;
    int itcount;
    int totitcount;
    int maxitcount;			
    int itzero;				
    int fcn;    

    /* output filename */
    char filename[STRMAX];			// Integration filename.
    int nout;                                   // number of output values.

    } solver_stat; 

	
#endif /*DEFH_*/
