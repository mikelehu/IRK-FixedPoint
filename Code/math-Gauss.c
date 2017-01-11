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
#include <GaussUserProblem.h> 
#include <GaussCommon.h>
#include <GaussCoefficients.h>
#include <time.h>
#include <stdbool.h>

/* Global variables */

int thread_count;


void mathGauss (int neq, int n, int ns, double t0, double tend,
                double *u0,long ulen, double *e0,long elen, 
                double h, double *rpar,long rlen, int *ipar,long ilen,
                int approximation, int threads,int algorithm,
                int rdigits1,int rdigits2,
                const char *myfilename,int sampling,int codfun)

{

/*------ declarations --------------------------------------------------------*/

    int i,is,initmean,totmean;
    solution u,u2; 
    solver_stat thestat,thestat2;
    gauss_method gsmethod,gsmethod2;
    ode_sys system;
    parameters params;
    toptions options,options2;
    int aux[2];
 
    u.uu = (val_type *)malloc(ulen*sizeof(val_type));
    u.ee = (val_type *)malloc(ulen*sizeof(val_type));
    u2.uu = (val_type *)malloc(ulen*sizeof(val_type));
    u2.ee = (val_type *)malloc(ulen*sizeof(val_type));

    options.rtol=malloc(ulen*sizeof(val_type));
    options.atol=malloc(ulen*sizeof(val_type));
    options2.rtol=malloc(ulen*sizeof(val_type));
    options2.atol=malloc(ulen*sizeof(val_type));

    clock_t clock0, clock1; 
    time_t  wtime0,wtime1;

    options.rdigits=rdigits1;    
    options2.rdigits=rdigits2;   


/* ----------- implementation  -----------------------------------------------*/     


    params.rpar =(val_type *)malloc(rlen*sizeof(val_type));
    params.ipar =(int *)malloc(ilen*sizeof(int));
  
    thread_count=threads;
    options.t0=t0;
    options.t1=tend;
    options.h=h;
    options2.h=h;
    options.algorithm=algorithm;
    options.sampling=sampling;
    options.approximation=approximation;
    gsmethod.ns=ns;
    gsmethod2.ns=gsmethod.ns;
    system.neq=neq;
    system.n=n;
    

    for (i=0; i<rlen; i++) params.rpar[i]=rpar[i];
    for (i=0; i<ilen; i++) params.ipar[i]=ipar[i];
    params.numrpar=rlen;


    select_odefun (codfun,&system);

    system.params.rpar=&params.rpar[0];
    system.params.ipar=&params.ipar[0];

    for (i=0; i<system.neq; i++)        
    {
          options.rtol[i]=RTOL;
          options.atol[i]=ATOL;
    }
  

    if (options.rdigits>0) options.mrdigits=pow(2,options.rdigits);
    if (options2.rdigits>0) options2.mrdigits=pow(2,options2.rdigits);

    if (options.algorithm==308 || options.algorithm==318)
    {
          GaussCoefficientsX(DIR_MATH,&gsmethod,&options);
          GaussCoefficientsX(DIR_MATH,&gsmethod2,&options2);
    }
    else 
    {    
          GaussCoefficients(DIR_MATH,&gsmethod,&options);
          GaussCoefficients(DIR_MATH,&gsmethod2,&options2);
    }


    strncpy(thestat.filename, myfilename,STRMAX); 
    InitStat(&system,&gsmethod,&thestat);

    for (i=0; i<neq;i++)
    { 
          u.uu[i]=u0[i];
          u.ee[i]=e0[i];					
    }


    wtime0= time(NULL);
    clock0= clock();


    select_gauss(&gsmethod, &gsmethod2, &u, &u2,
           &system, &options, &options2,
           &thestat, &thestat2);

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

    switch (options.algorithm)
    {  
    case 1: case 2:  case 11:  case 12:  case 301: case 308: case 309:
    case 311: case 319: case 303: 

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

    break;



    case 321: case 322: 

          aux[0]=thestat.totitcount;
          aux[1]=thestat2.totitcount;

          if (!MLPutFunction(stdlink, "List",14)) MLErrorMessage(stdlink);
          if (!MLPutReal64( stdlink,thestat.MaxDE)) MLErrorMessage(stdlink); 
          if (!MLPutReal64List( stdlink, u.uu, neq)) MLErrorMessage(stdlink); 
          if (!MLPutReal64List( stdlink, u.ee, neq))  MLErrorMessage(stdlink);
          if (!!MLPutInteger32( stdlink,thestat.stepcount))
              MLErrorMessage(stdlink);
          if (!!MLPutInteger32List (stdlink,aux,2)) MLErrorMessage(stdlink); 
          if (!!MLPutInteger32( stdlink,thestat.maxitcount))
              MLErrorMessage(stdlink); 
          if (!!MLPutInteger32( stdlink,thestat.itzero))
              MLErrorMessage(stdlink);
          if (!!MLPutInteger32( stdlink,thestat.fcn)) MLErrorMessage(stdlink); 
          if (!MLPutReal64( stdlink,(float) (clock1 - clock0)/CLOCKS_PER_SEC)) 
              MLErrorMessage(stdlink); 
          if (!MLPutReal64( stdlink,(float) (wtime1 - wtime0)))
              MLErrorMessage(stdlink); 
          if (!MLPutInteger32( stdlink,thestat.convergence))
              MLErrorMessage(stdlink);   
          if (!MLPutInteger32( stdlink,thestat.nout)) MLErrorMessage(stdlink); 
          if (!MLPutInteger32( stdlink,totmean)) 
             MLErrorMessage(stdlink); 
          if (!MLPutInteger32( stdlink,thestat.totitcountzero)) 
             MLErrorMessage(stdlink);            
          if(! MLEndPacket(stdlink)) MLErrorMessage(stdlink);
          if(! MLFlush(stdlink)) MLErrorMessage(stdlink);
    break;

    default:
           printf("Error: incorrect algorithm\n");
    break;

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

}


int main(int argc, char *argv[])
{
     return MLMain(argc, argv);
}


