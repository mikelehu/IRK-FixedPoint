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

   def.c:               Parameters and general definitions we use in the code.
                        You must specify math-functions of your Odefun. 

   Terminal-IRKFP.c:    An example to show how to run an integration from the Unix terminal.
   math-IRKFP.c:	An example to show how to run an integration from mathematica. 

   GaussCoefficients.c: Function to load Butcher tableau of the Implicit Runge-Kutta method and 
                        the coefficients of the  reformulation as explained in the article.
                        We offer Gauss method coefficients for s=6, s=8, s=16.

   Common-IRKFP.c:      Numerical integration method. We will find these functions:

	IRKFP () : 	         IRK Fixed-Point integration (main function).
	Fixed_point_Step ():     Fixed iteration Step.	
	General_FP_It():         General fixed point iteration.	
	Partitioned_FP_It():     Partitioned fixed point iteration.		
        Default_Stage_init:      Initialization of Stages, Y_{n,i}=y_n	
        Interpolated_Stage_init	 Initialization of stages, Y_{n,i}=G(Y_{n-1,i})		
											
	MyOutput():	         An output function defined by the user.

	NormalizedDistance():    Check convergence of the iteration.
	StopCriterion ():        Stopping criterion.
	RemoveDigits():          Rounding a floating point number with p-r significant binary digits.
        CompensatedSummation():  Function to compute y_{n+1}=y_{n}+ sum_{i} L_{n,i}				
							
   Problems.c:  Examples of ODE systems and Hamiltonians: Double pendulum and N-Body problems.
														
	Ode1()=OdePendulum():	
	Ham1()=HamPendulum():	
	Ode2()=OdeNBody():
	Ham2()=HamNBody():

   math-IRKFP.tm:		 Mathlink file.
			

*********************************************************************************
ARGUMENTS:

   You will have to specify next options for the numerical integration:
     IRKFP (t0,t1,h, &gsmethod, &u, &system, &options, &thestat);
      
       t0,t1: interval of numerical integration.
       h:     stepsize.

       gsmethod->ns: number of stages of Implicit Runge-Kutta Gauss collocation method.
       gsmethod->m, hc, hb: coefficients that define Implicit Runge-Kutta method.
       hsmethod->nu: interpolation coefficients to initialize stages.

       u->uu[neq], u->ee[neq]: initial values for the numerical integration.

       system->neq: dimension of the differential equation.
       system->eda: name of the function that compute differential equation.
       system->ham: name of the Hamiltonian.
       system->rpar[]: real paramaters of the ODE system.
       system->ipar[]: integer parameters of the ODE system. 

       options.rtol[neq]:   relative tolerance.
       options.atol[neq]:   absolute tolerance.
       options.TheOutput:   name of the function defined by the user 
                         that execute at each step of the integration.
       options.filename:    name of the filename to write the 
                         numerical solution of the integration.
       options.iteration:   kind of numerical iteration, general or partitioned.
       options.StageInitFn: name of the functon for the initialization of stages.
       options.rdigits:    number of binary digits we remove from L_{n,i}
       options.sampling:    we sample the numerical results once every "sampling" steps.
       options.filename:    name of the output binary filename. 

   thread_count: number of threads for parallel computation. 


*********************************************************************************
PARAMETERS (file: def.h):

   You can specify next parameters:

   PARALLEL: specify to run a parallel execution.

   MAXIT:  maximum number of fixed point iterations 
   RTOL,ATOL: fixed point iteration default tolerances. 


   #define DIR_TERM :  // Path Coefficients for terminal executions  
   #define DIR_MATH :  // Path Coefficients for mathematica executions

*********************************************************************************
RESULTS:

   options.TheOutput: user defined function.
   options.filename:  output binary file.

********************************************************************************
INSTALATION (Ubuntu 16.04 lts):

   Two ways to execute :

     Unix terminal execution:
          make term-IRKFP
          ./Terminal-IRKFP.exe

      Mathematica (double)
          make math-IRKFP
          Execution from mathematica notebook (see examples).


********************************************************************************
EXAMPLES:

    NCDP: Mathematica notebook integration of Double Pendulum problem (non-chaotic).
    CDP : Mathematica notebook integration of Double Pendulum problem (chaotic).
    OSS : Mathematica notebook integration of Outer solar system problem.

*********************************************************************************

    
    



