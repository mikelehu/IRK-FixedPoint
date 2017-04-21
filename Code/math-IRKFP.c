/*----------------------------------------------------------------------------*/
/*									      */
/*                                math-Gauss.c				      */
/*									      */
/* ---------------------------------------------------------------------------*/

#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <mathlink.h>
#include <def.h>
#include <stdio.h>
#include <Problems.h> 
#include <Common-IRKFP.h>
#include <GaussCoefficients.h>
#include <time.h>
#include <stdbool.h>

/* Global variables */

int thread_count;


void mathIRKFP (int neq, int ns, double t0, double t1,
                double *u0,long ulen, double *e0,long elen, 
                double h, double *rpar,long rlen, int *ipar,long ilen,
                int approximation, int threads,int iteration,
                int rdigits,
                const char *myfilename,int sampling,int codfun)

{

/*------ declarations --------------------------------------------------------*/

    int i,is,initmean,totmean;

    gauss_method gsmethod;
    solution u; 
    toptions options;
    ode_sys system;
    solver_stat thestat;
   
    clock_t clock0, clock1; 
    time_t  wtime0,wtime1;         


/* ----------- implementation  --------------------------------------------*/     

    gsmethod.ns=ns;

/* ----------- problem parameters-----------------------------------------*/

    system.neq=neq;

    switch (codfun)
    { 
     case 1: 
           system.f = Ode1;
           system.ham= Ham1;                   
     break;

     case 2: 
           system.f = Ode2;
           system.ham= Ham2;
     break;

     default:
           printf("error. codfun\n");
     break;
     
    }

/* Parameters of the system */

    system.params.rpar =(val_type *)malloc(rlen*sizeof(val_type));
    system.params.ipar =(int *)malloc(ilen*sizeof(int));
    system.params.numrpar=rlen;
    system.params.numipar=ilen;
    for (i=0; i<rlen; i++) system.params.rpar[i]=rpar[i];
    for (i=0; i<ilen; i++) system.params.ipar[i]=ipar[i];

/* ----------- Initial values--------------------------------------------*/

    u.uu = (val_type *)malloc(ulen*sizeof(val_type));
    u.ee = (val_type *)malloc(ulen*sizeof(val_type));

    for (i=0; i<neq;i++)
    { 
          u.uu[i]=u0[i];
          u.ee[i]=e0[i];					
    }


/* ----------- Integration options---------------------------------------*/


    options.rtol=malloc(ulen*sizeof(val_type));
    options.atol=malloc(ulen*sizeof(val_type));

    for (i=0; i<system.neq; i++)        
    {
          options.rtol[i]=RTOL;
          options.atol[i]=ATOL;
    }


    options.TheOutput=MyOutput;
    strncpy(options.filename, myfilename,STRMAX); 

    switch (iteration)
    {
     case 0:
       options.iteration=General_FP_It;    
     break;

     case 1:
       options.iteration=Partitioned_FP_It;
     break;

     default:
       printf("error. iteration\n");
     break;
     
    }   

   
    switch (approximation)
    { 
     case 0: 
           options.StageInitFn=Default_Stage_init;                  
     break;

     case 1: 
           options.StageInitFn=Interpolated_Stage_init;
     break;

     default:
           printf("error. Stage initialization\n");
     break;
     
    }    

    options.rdigits=rdigits; 
    options.mrdigits=pow(2,options.rdigits);

/* ----------- Integration interval----------------------------------------*/
  
    options.sampling=sampling;

    
/* ----------- Integration of the problem  ------------------------------------*/

    thread_count=threads;      

    GaussCoefficients(DIR_MATH,&gsmethod,h);

    InitStat(&system,&gsmethod,&thestat);


    wtime0= time(NULL);
    clock0= clock();


    IRKFP (t0,t1,h, &gsmethod, &u, &system, &options, &thestat);  

    clock1=clock();
    wtime1= time(NULL);

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

    if (!MLPutFunction(stdlink, "List",14)) MLErrorMessage(stdlink);
    if (!MLPutReal64( stdlink,thestat.MaxDE)) MLErrorMessage(stdlink); 
    if (!MLPutReal64List( stdlink, u.uu, neq))  MLErrorMessage(stdlink); 
    if (!MLPutReal64List( stdlink, u.ee, neq))  MLErrorMessage(stdlink); 
    if (!!MLPutInteger32( stdlink,thestat.stepcount))
        MLErrorMessage(stdlink);
    if (!!MLPutInteger32( stdlink,thestat.totitcount))
        MLErrorMessage(stdlink); 
    if (!!MLPutInteger32( stdlink,thestat.maxitcount))
        MLErrorMessage(stdlink); 
    if (!!MLPutInteger32( stdlink,thestat.itzero))
        MLErrorMessage(stdlink); 
    if (!!MLPutInteger32( stdlink,thestat.fcn))
        MLErrorMessage(stdlink); 
    if (!MLPutReal64( stdlink,(float) (clock1 - clock0)/CLOCKS_PER_SEC))
        MLErrorMessage(stdlink);
    if (!MLPutReal64( stdlink,(float) (wtime1 - wtime0)))
        MLErrorMessage(stdlink); 
    if (!MLPutInteger32( stdlink,thestat.convergence))
        MLErrorMessage(stdlink); 
    if (!MLPutInteger32( stdlink,thestat.nout))
        MLErrorMessage(stdlink);  
    if (!MLPutInteger32( stdlink,totmean)) 
        MLErrorMessage(stdlink); 
    if (!MLPutInteger32( stdlink,thestat.totitcountzero)) 
        MLErrorMessage(stdlink);    
    if(! MLEndPacket(stdlink)) MLErrorMessage(stdlink);
    if(! MLFlush(stdlink))  MLErrorMessage(stdlink); 

   
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

}


int main(int argc, char *argv[])
{
     return MLMain(argc, argv);
}


