****************************************************************************
Readme.txt
****************************************************************************
     Gauss Implicit Runge-Kutta implementation. 
     Fixed-Point iteration.

     version: 1 (04-02-2016).

     Article: "Reducing and monitoring round-off error propagation for sympletic
              implicit Runge-Kutta schemes" (2016). 

********************************************************************************
CONTENTS:

   prec.c:		Precision.
   def.c:               Parameters and general definitions we use in the code.
                        You must specify math-functions of you Odefun. 

   GaussTerminal.c:     An example to show how to call the numerical integration.
   math-Gauss.c:	An auxiliar file to call from mathematica (double). 
   math-GaussPar.c:     An auxiliar file to call from mathematica (two integration in parallel).
   quad-Gauss.c:	An auxiliar file to call from mathematica (quadruple).

   GaussInitData.c:     Contain "InitialData" with the initial values for double pendulum problem.

   GaussCoefficients.c: Coefficients mij,bi,ci that define Inplicit Runge-Kutta method.

   GaussCommon.c:       Numerical integration method. We will find these functions:

	RKG () : 	   IRK method.
	RKG2():		   Both primary and secondary sequential for round-off error estimation.
	Fixed_point_it (): Fixed iteration method.	
	It_Jacobi():								
	It_Seidel():	
	Yi_init();												
	TheOutput():	   Output function for RKG ().
	TheOutput2():	   Output function for RKG-2 ().

	NormalizedDistance():
	UpdateDMin ():
	RemoveDigitsFcn():				
							
   GaussUserProblem.c:  Double pendulum and N-Body ode system:
														
	OdePendulum():	
	HamPendulum():	
	OdeNBody():
	HamNBody():

   math-Gauss.tm:	Mathlink file.
   quad-Gauss.tm:	Mathlink file.
   math-GaussPar.tm:	Mathlink file.				

*********************************************************************************
OPTIONS:

   You will have to specify next options for the numerical integration:

   ns= number of stages of Inplicit Runge-Kutta method.
   eda= name of differential equation.
   algorith= You can execute four differents kind of itegrations:
	=  1 Jacobi fixed point iteration method. 
	= 11 Seidel fixed point iteration method.
	= 21 Both integrations execute sequentially (Jacobi).
        = 22 Both integrations execute sequentially (Seidel).
   
   h= stepsize.
   sampling: for each span of day we need the output.
   filename = output binary filename. 

   thread_count: number of thread out computer have. ???

   note: neq, t0, tf and some others options of the problems are 
         initialized with initial values.

*********************************************************************************
PARAMETERS:

   You can specify next parameters:

   IOUT: default form (when we need intermedi values). 
   PARALLEL: we want parallel execution.

   MAXIT 50 :  maximum number of fixed point iterations 
   MAXKSW 10:  max number of steps to change second integration initialization mode.
   RTOL,ATOL: fixed point iteration tolerance. 
   PREC:
       1=DOUBLEPRECISION
       2=QUADRUPLEPRECISION
       3=FLOAT.

*********************************************************************************
RESULTS:

   filename (output binary format).

********************************************************************************
INSTALATION:
   Three options:

     Terminal execution:
          make term-Gauss
          ./GaussTerminal.exe

      Mathematica (double)
          make math-Gauss
          Execution from mathematica notebook.

      Mathematica (quadruple)
          make quad-Gauss
          Execution from mathematica notebook.

      Mathematica (two integration in parallel for estimating error).
          make math-GaussPar
          Execution from mathematica notebook. 

********************************************************************************
EXAMPLES:
    RKGC_DoublePendulum-I.nb: Mathematica notebook integration of Double Pendulum problem (non-chaotic).
    RKGC_DoublePendulum-II.nb: Mathematica notebook integration of Double Pendulum problem (chaotic).
    RKGC_NBody.nb: Mathematica notebook integration of N9-Body problem.

*********************************************************************************
QUADRUPLE PRECISION:
    We have used libquadmath GCC Quad-Precision Math Library.
    https://gcc.gnu.org/onlinedocs/libquadmath/.
    
    



