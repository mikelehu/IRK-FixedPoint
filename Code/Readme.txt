****************************************************************************
Readme.txt
****************************************************************************
     Gauss Implicit Runge-Kutta implementation. 
     Fixed-Point iteration.

     version: 1.1 (21-03-2017).

     Article: "Reducing and monitoring round-off error propagation for sympletic
               implicit Runge-Kutta schemes" (2017). 

********************************************************************************
CONTENTS:

   prec.c:		Precision.
   def.c:               Parameters and general definitions we use in the code.
                        You must specify math-functions of your Odefun. 

   GaussTerminal.c:     An example to show how to call the numerical integration.
   math-Gauss.c:	An auxiliar file to call from mathematica (double precision). 

   GaussInitData.c:     Contain "InitialData" with the initial values for some problems.

   GaussCoefficients.c: Coefficients mij,bi,ci that define Inplicit Runge-Kutta method.

   GaussCommon.c:       Numerical integration method. We will find these functions:

	RKG () : 	   IRK method.
	RKG2():		   Both primary and secondary sequential for round-off error estimation.
	Fixed_point_it (): Fixed iteration method.	
	It_Jacobi():	   Standard fixed point iteration.							
	Yi_init();	   Initialization functions for Yi stages.
											
	TheOutput():	   Output function for RKG ().
	TheOutput2():	   Output function for RKG-2 ().

	NormalizedDistance(): check of convergence of the iteration.
	StopCriterion ():   : Stopping criterion.
	RemoveDigitsFcn()   : Rounding a floating point number with p-r significant binary digits.				
							
   GaussUserProblem.c:  Double pendulum and N-Body ode system:
														
	Ode1()=OdePendulum():	
	Ham1()=HamPendulum():	
	Ode2()=OdeNBody():
	Ham2()=HamNBody():

   math-Gauss.tm:	Mathlink file.
			

*********************************************************************************
OPTIONS:

   You will have to specify next options for the numerical integration:

   ns= number of stages of Inplicit Runge-Kutta method.
   eda= name of differential equation.

   algorithm= You must specify one of the next implementations options:
	=  1 Standard fixed point iteration method. 
	=  2 Both integrations execute sequentially.

   
   h= stepsize.
   sampling: we sample the numerical results once every "sampling" steps.
   filename = output binary filename. 

   thread_count: number of threads for parallel computation. 

   note: neq, t0, tf and some others options of the problems are 
         initialized with initial values.

*********************************************************************************
PARAMETERS (file: def.h):

   You can specify next parameters:

   IOUT: default form (enable otuput binary file). 
   PARALLEL: specify to run a parallel execution.

   MAXIT 50 :  maximum number of fixed point iterations 
   MAXKSW 10:  max number of steps to change second integration initialization mode.
   RTOL,ATOL: fixed point iteration tolerances. 
   PREC:
       1=DOUBLEPRECISION

   #define DIR_TERM :  // Path Coefficients for terminal executions  
   #define DIR_MATH :  // Path Coefficients for mathematica executions

*********************************************************************************
RESULTS:

   filename (output binary format).

********************************************************************************
INSTALATION (Ubuntu 16.04 lts):

   Two ways to execute :

     Terminal execution:
          make term-Gauss
          ./GaussTerminal.exe

      Mathematica (double)
          make math-Gauss
          Execution from mathematica notebook (see examples).


********************************************************************************
EXAMPLES:

    NCDP: Mathematica notebook integration of Double Pendulum problem (non-chaotic).
    CDP : Mathematica notebook integration of Double Pendulum problem (chaotic).
    OSS : Mathematica notebook integration of Outer solar system problem.

*********************************************************************************

    
    



