/*------------------------------------------------------------------------------*/
/*										*/
/*                                GaussTerminal.c				*/
/*										*/
/* -----------------------------------------------------------------------------*/

#include <GaussTerminal.h>

/* Global variables */

int thread_count=1;

int main()
{

/*------ declarations --------------------------------------------------*/    

     int i;   
	
     gauss_method gsmethod,gsmethod2;
     solution u,u2;
     toptions options,options2;
     parameters params;
     ode_sys system;   
     solver_stat thestat,thestat2;

     clock_t clock0, clock1; 
     time_t  wtime0,wtime1;

     params.rpar =(val_type *)malloc(MAXPARAM*sizeof(val_type));
     params.ipar =(int *)malloc(MAXPARAM*sizeof(int));

     u.uu = (val_type *)malloc(MAXNEQ*sizeof(val_type));
     u.ee = (val_type *)malloc(MAXNEQ*sizeof(val_type));
     u2.uu = (val_type *)malloc(MAXNEQ*sizeof(val_type));
     u2.ee = (val_type *)malloc(MAXNEQ*sizeof(val_type));

     options.rtol=malloc(MAXNEQ*sizeof(val_type));
     options.atol=malloc(MAXNEQ*sizeof(val_type));
  
#    if PREC ==2  //QUADRUPLEPRECISION
     int n;
     int width = 46;
     char buf[128];
#    endif

/* ----------- implementation  ---------------------------------------*/     


/* ----------- integration parameters---------------------------------*/  

     gsmethod.ns = 6;     				 //	 Stages.
     options.h = POW(2,-7);	       			 //	 Stepsize. 
     options.sampling=1;   			 	 //      

     system.f = Ode1;					 //	 Odefun (GaussUserProblem.c: OdePendulum,OdeNBody)
     system.ham= Ham1;					 //      Hamiltonian (GaussUserProblem.c: HamPendulum,HamNBody)

     system.problem =1; 				 //	 Initial values (GaussInitData.c).

     options.approximation=1;   			 //      Approximation: Y^[0] (GaussCommon.c/Yi_init()).
	
     strncpy(thestat.filename, "AllSolve1.bin",STRMAX);   //     Output filename.

     options.algorithm=1;    				 //	 1=Jacobi;  11=Seidel;    
     options.rdigits=0;       
     options2.rdigits=3;

/* ----------- execution  ------------------------------------------*/

     printf("Begin execution \n");
     printf("method=%i, problem=%i, algorithm=%i\n",gsmethod.ns,system.problem,options.algorithm);
     printf("approximation=%i,sampling=%i\n",options.approximation,options.sampling);

#if PREC ==2  //QUADRUPLEPRECISION
     printf("options.h=");
     n = quadmath_snprintf(buf, sizeof buf, "%+-#*.30Qe", width, options.h);
     if ((size_t) n < sizeof buf) printf("%s\n",buf);
#else         //DOUBLEPRECISION
     printf("options.h=%lg\n", options.h);        
#endif 
     printf("----------------\n");
        
     InitialData (&options,&u,&params,&system);                 
     system.params.rpar=&params.rpar[0];
     system.params.ipar=&params.ipar[0];
     
     for (i=0; i<system.neq; i++)
     {
          options.rtol[i]=RTOL;
          options.atol[i]=ATOL;
     }

     if (options.rdigits>0) options.mrdigits=pow(2,options.rdigits);
     if (options2.rdigits>0) options2.mrdigits=pow(2,options2.rdigits);

     GaussCoefficients(&gsmethod,&options);
   
     InitStat(&system,&gsmethod,&thestat);	
   

     print_u(system.neq, u.uu);

     wtime0= time(NULL);
     clock0= clock();

     select_gauss(&gsmethod, &gsmethod2, &u, &u2,&system,
                  &options, &options2, &thestat, &thestat2);
  
     clock1=clock();
     wtime1= time(NULL);

     printf("End execution \n");
        
     if (thestat.convergence==SUCCESS)
     {
           printf("Execution Correct\n");
           printf("convergence=%i.  (=0 correct;=-1 incorrect)\n",thestat.convergence);        

           print_u (system.neq,u.uu);
#if PREC ==2  //QUADRUPLEPRECISION
           printf("Max-DE");
           n = quadmath_snprintf(buf, sizeof buf, "%+-#*.30Qe", width, thestat.MaxDE);
           if ((size_t) n < sizeof buf) printf("%s\n",buf);
#else  // DOUBLEPRECISION    
           printf("Energy MaxDE=%.20lg\n",thestat.MaxDE);
#endif                      

           printf ("\nCPU time:%lg\n", (double) (clock1 - clock0)/CLOCKS_PER_SEC);
	   printf ("Elapsed wall clock time: %ld\n", (wtime1 - wtime0));
           printf ("Elapsed wall clock time: %lg\n", difftime(wtime1, wtime0));

           printf("\n");
           
           printf("stepcount=%i\n",thestat.stepcount);
           printf("nout=%i\n",thestat.nout);
           
           printf("fixed-point iterations=%i\n",thestat.totitcount);
           printf("max fixed-point iterations=%i\n",thestat.maxitcount);        
           printf("deltaZero iterations=%i\n",thestat.itzero); 
 
           printf("\n");

           printf ("function ebaluations\n");
           printf("fcn=%i\n",thestat.fcn);
           printf ("\n");

     }
     else
     {
           printf("Execution InCorrect\n");
           printf("convergence=%i.  (=0 correct;=-1 incorrect)\n",thestat.convergence);
     }
 
     free(params.rpar); free(params.ipar);
     free(u.uu); free(u.ee);
     free(u2.uu); free(u2.ee);
    
     exit(0);

}





