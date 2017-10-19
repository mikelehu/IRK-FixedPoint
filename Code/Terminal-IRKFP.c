/* ---------------------------------------------------------------------------*/
/*								              */
/*                                GaussTerminal.c			      */
/*								              */
/* ---------------------------------------------------------------------------*/

#include <Terminal-IRKFP.h>

/* Global variables */

int thread_count=1;

int main()
{

/*------ declarations --------------------------------------------------------*/    

     int i,is,initmean,totmean;   
	
     gauss_method gsmethod;
     solution u;
     toptions options;
     ode_sys system;   
     solver_stat thestat;
     val_type h, t0, t1;

     clock_t clock0, clock1; 
     time_t  wtime0,wtime1;

     /* Auxiliar variables for the problem */
     val_type m1,m2,l1,l2,g;
     val_type C1,C2,C3,C4,C5,C6,C7;


/* ----------- implementation  --------------------------------------------*/     


/* ----------- integration parameters-------------------------------------*/  

     gsmethod.ns = 6;     	 // Stages.

/* ----------- problem parameters-----------------------------------------*/

/* Double pendulum  */ 
	
     system.neq=4;		 	
     system.f = Ode1;		// User defined ODE system: \dot{u}=f(t,u)
                                // The function call must be: f(neq, t, u, fu, params)
     system.ham= Ham1;          // User defined Hamiltonian of the system
                                // The function call must be: H(neq,u, params)     
     system.cod[0]=0;

/* Parameters of the system */

     m1=1.;
     m2=1.;
     l1=1.;
     l2=1.;
     g=9.8;

     C1=-l1*l1*(m1+m2);
     C2=-l2*l2*m2;
     C3=-2*l1*l2*m2;
     C4=-2*l1*l1*l2*l2*m2*m1-l1*l1*l2*l2*m2*m2;
     C5=-l1*l1*l2*l2*m2*m2;
     C6=g*l1*(m1+m2);
     C7=g*l2*m2;

     system.params.rpar =(val_type *)malloc(7*sizeof(val_type));
     system.params.numrpar=7;
     system.params.ipar =(int *) 0;
     system.params.numipar=0;    
     system.params.rpar[0]=C1;
     system.params.rpar[1]=C2;	
     system.params.rpar[2]=C3;
     system.params.rpar[3]=C4;
     system.params.rpar[4]=C5;
     system.params.rpar[5]=C6;
     system.params.rpar[6]=C7;

/* ----------- Initial values-----------------------------------------------*/

     u.uu = (val_type *)malloc(system.neq*sizeof(val_type));
     u.ee = (val_type *)malloc(system.neq*sizeof(val_type));

    //u0=(q,p) 
     u.uu[0]=11./10;       // q0
     u.uu[1]=-11./10; 	   // q1
     u.uu[2]=13873./5000;  // p0
     u.uu[3]=13873./5000;  // p1

    //e0
     u.ee[0]=-8.8817841970012523234e-17;
     u.ee[1]= 8.8817841970012523234e-17;
     u.ee[2]=4.4764192352886311710e-17; 
     u.ee[3]=4.4764192352886311710e-17;      


/* ----------- Integration options------------------------------------------*/

    options.rtol=malloc(system.neq*sizeof(val_type));
    options.atol=malloc(system.neq*sizeof(val_type));

    for (i=0; i<system.neq; i++)
     {
          options.rtol[i]=RTOL;
          options.atol[i]=ATOL;
     }
    
     options.TheOutput=MyOutput;  // Function to be called after each step.
                                  // Defined by the user.
     strncpy(options.filename, "Output.bin",STRMAX);   // Output filename.

     options.iteration=General_FP_It;         // Partitioned_FP_It
     options.StageInitFn=Default_Stage_init;
     options.rdigits=0;         	      // =0 Double Precision.
     options.mrdigits=pow(2,options.rdigits); // If rdigits=0 it is not used.

/* ----------- Integration interval--------------------------------------------*/

     h = POW(2,-7);
     t0=0.;
     t1=pow(2,12);  
     options.sampling=POW(2,10);

/* ----------- Integration of the problem  ------------------------------------*/

     printf("Begin integration \n"); 

     GaussCoefficients(DIR_TERM,&gsmethod,h);    // Initialize Coefficients of the method

     InitStat(&system,&gsmethod,&thestat);       // Initialize the stat

     wtime0= time(NULL);
     clock0= clock();


     IRKFP (t0,t1,h, &gsmethod, &u, &system, &options, &thestat);  

  
     clock1=clock();
     wtime1= time(NULL);

     printf("End integration \n"); 

/* --------- Results ---------------------------------------------------------*/ 


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
 
    free(system.params.rpar); 
    free(u.uu);
    free(u.ee);

    free(options.rtol); 
    free(options.atol);


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





