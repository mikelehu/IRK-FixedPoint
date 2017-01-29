/* ---------------------------------------------------------------------------*/
/*								              */
/*                                GaussTerminal.c			      */
/*								              */
/* ---------------------------------------------------------------------------*/

#include <GaussTerminal.h>

/* Global variables */

int thread_count=1;

int main()
{

/*------ declarations --------------------------------------------------------*/    

     int i,is,initmean,totmean,codfun;   
	
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
     options2.rtol=malloc(MAXNEQ*sizeof(val_type));
     options2.atol=malloc(MAXNEQ*sizeof(val_type));
  

/* ----------- implementation  -----------------------------------------------*/     


/* ----------- integration parameters-----------------------------------------*/  

     gsmethod.ns = 6;     	 // Stages.
     gsmethod2.ns=gsmethod.ns;

     strncpy(thestat.filename, "Output.bin",STRMAX);   // Output filename.

     options.approximation=0;  	 // Approximation: Y^[0] (GaussCommon.c/Yi_init()).

     options.algorithm=1;    	 // (1=Fixed point iteration,2=Error estimation).    
     options.rdigits=0;       
     options2.rdigits=3;

/* ----------- problem parameters-----------------------------------------*/

/* Double pendulum : NCDP or CDP */ 
			 	
     codfun = 1;                 // Odefun
     system.problem =1; 	 // NCDP - Initial values (GaussInitData.c).
     system.problem =2;          // CDP  - Initial values (GaussInitData.c).


/* Outer solar system            */

     codfun = 2;                 // Odefun
     system.problem =3;          // OSS  - Initial values (GaussInitData.c).


/* ----------- execution  ----------------------------------------------------*/

     printf("Begin execution \n");
        
     InitialData (&options,&u,&params,&system);

     select_odefun (codfun,&system);
     options2.h = options.h;
                 
     system.params.rpar=&params.rpar[0];
     system.params.ipar=&params.ipar[0];

     printf("method=%i, problem=%i, algorithm=%i\n",
             gsmethod.ns,system.problem,options.algorithm);
     printf("approximation=%i,sampling=%i\n",
             options.approximation,options.sampling);

     printf("options.h=%lg\n", options.h);        
     printf("----------------\n");
     
     for (i=0; i<system.neq; i++)
     {
          options.rtol[i]=RTOL;
          options.atol[i]=ATOL;
     }


     if (options.rdigits>0) options.mrdigits=pow(2,options.rdigits);
     if (options2.rdigits>0) options2.mrdigits=pow(2,options2.rdigits);

     GaussCoefficients(DIR_TERM,&gsmethod,&options); 
     GaussCoefficients(DIR_TERM,&gsmethod2,&options2); 

     InitStat(&system,&gsmethod,&thestat);

     wtime0= time(NULL);
     clock0= clock();

     select_gauss(&gsmethod, &gsmethod2, &u, &u2, &system,
                  &options, &options2, &thestat, &thestat2);
  
     clock1=clock();
     wtime1= time(NULL);

/* --------- Results ---------------------------------------------------------*/
     printf("End execution \n");   


    totmean=100;
    if (thestat.stepcount>1)
       for (is=0; is<gsmethod.ns; is++)
       {
            for (i=0; i<system.neq; i++) 
            {
                  initmean=(int)thestat.initqlty[is*system.neq+i]/(thestat.stepcount-1);
                  if (initmean < totmean) totmean=initmean;
            }               
       }
     else totmean=0;



     if (thestat.convergence==SUCCESS)
     {
           printf("Execution Correct\n");
           printf("convergence=%i.  (=0 correct;=-1 incorrect)\n",
                   thestat.convergence);        

           print_u (system.neq,u.uu); 
           printf("Energy MaxDE=%.20lg\n",thestat.MaxDE);                     

           printf ("\nCPU time:%lg\n", 
                  (double) (clock1 - clock0)/CLOCKS_PER_SEC);
	   printf ("Elapsed wall clock time: %ld\n", (wtime1 - wtime0));
           printf ("Elapsed wall clock time: %lg\n", difftime(wtime1, wtime0));

           printf("\n");
           
           printf("stepcount=%i\n",thestat.stepcount);
           printf("nout=%i\n",thestat.nout);
           
           printf("fixed-point iterations=%li\n",thestat.totitcount);
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
           printf("convergence=%i.  (=0 correct;=-1 incorrect)\n",
                   thestat.convergence);
     }
 
    free(params.rpar); free(params.ipar);
    free(u.uu); free(u.ee);
    free(u2.uu); free(u2.ee);

    free(options.rtol); free(options.atol);
    free(options2.rtol); free(options2.atol);

    free(gsmethod.m); 
    free(gsmethod.a);
    free(gsmethod.b); 
    free(gsmethod.hb);      
    free(gsmethod.c); 
    free(gsmethod.hc); 
    free(gsmethod.nu); 
    free(gsmethod.orderedindices);

    free(thestat.z); 
    free(thestat.li);
    free(thestat.fz);
    free(thestat.initqlty);
    free(thestat.zit0);
    
    exit(0);

}





