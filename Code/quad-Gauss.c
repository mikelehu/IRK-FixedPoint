/*------------------------------------------------------------------------------*/
/*										*/
/*                                quad-Gauss.c					*/
/*										*/
/* -----------------------------------------------------------------------------*/

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

void quadGauss(int neq, int n, int ns, double t0, double tend,
	const unsigned char *quadak, int qluz, const unsigned char *elsak, int elluz,
	double h, 
	const unsigned char *realpars, int numrealpars, int *ipar, long ilen, 
	int approximation, int threads, int algorithm, 
        int rdigits1,int rdigits2, 
	const char *myfilename, int sampling, int codfun)
{
//------ declarations --------------------------------------------------

    int i;
    solution u,u2; 
    solver_stat thestat,thestat2;
    gauss_method gsmethod,gsmethod2;
    ode_sys system;
    parameters params;
    toptions options,options2;
    int aux[2];
 

    u.uu = (__float128 *)malloc(qluz);
    u.ee = (__float128 *)malloc(elluz);
    u2.uu = (__float128 *)malloc(qluz);
    u2.ee = (__float128 *)malloc(elluz);

    options.rtol=(__float128 *)malloc(qluz*sizeof(__float128));
    options.atol=(__float128 *)malloc(qluz*sizeof(__float128));
    options2.rtol=(__float128 *)malloc(qluz*sizeof(__float128));
    options2.atol=(__float128 *)malloc(qluz*sizeof(__float128));

    clock_t clock0, clock1; 
    time_t  wtime0,wtime1;

    options.rdigits=rdigits1; 
    options2.rdigits=rdigits2;      


// ----------- implementation  --------------------------------------    


    params.rpar =(__float128 *)malloc(numrealpars);
    params.ipar =(int *)malloc(ilen*sizeof(int));
  
    thread_count=threads;
    options.t0=t0;
    options.t1=tend;
    options.h=h;
    options.algorithm=algorithm;
    options.sampling=sampling;
    options.approximation=approximation;
    gsmethod.ns=ns;
    system.neq=neq;
    system.n=n;

    memcpy((void *)params.rpar,(const void *)realpars,(size_t)numrealpars);
    for (i=0; i<ilen; i++) params.ipar[i]=ipar[i];

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

    case 3: 
          system.f = Ode3;
          system.ham= Ham3;               
    break;

    case 4: 
          system.f = Ode4;
          system.ham= Ham4;                  
    break;

    case 5: 
          system.f = Ode5;
          system.ham= Ham5;               
    break;

    case 6: 
          system.f = Ode6;
          system.ham= Ham6;               
    break;

    case 7: 
          system.f = Ode7;
          system.ham= Ham7;                   
    break;

    case 8: 
          system.f = Ode8;
          system.ham= Ham8;                  
    break;

    case 9: 
          system.f = Ode9;
          system.ham= Ham9;                 
    break;

    case 10: 
          system.f = Ode10;
          system.ham= Ham10;                   
    break;

    case 11: 
          system.f = Ode11;
          system.ham= Ham1;
    break;

    case 12: 
          system.f = Ode12;
          system.ham= Ham2;
    break;

    default:
          printf("error. codfun\n");
    break;
    }
 
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
    strncpy(thestat.filename, myfilename,STRMAX); 
    InitStat(&system,&gsmethod,&thestat);

    memcpy((void *)u.uu,(const void *)quadak,(size_t)qluz);
    memcpy((void *)u.ee,(const void *)elsak,(size_t)elluz);
  
    wtime0= time(NULL);
    clock0= clock();

    select_gauss(&gsmethod, &gsmethod2, &u, &u2,
                 &system, &options, &options2,
                 &thestat, &thestat2);

    clock1=clock();
    wtime1= time(NULL);

    switch (options.algorithm)
    {  
    case 1: case 11: case 9: case 19:

          if (!MLPutFunction(stdlink, "List",12)) MLErrorMessage(stdlink);
          if (!MLPutReal64( stdlink,thestat.MaxDE)) MLErrorMessage(stdlink); 
	  if (!MLPutByteString( stdlink, (const unsigned char *)u.uu,qluz)) MLErrorMessage(stdlink);
	  if (!MLPutByteString( stdlink, (const unsigned char *)u.ee,elluz)) MLErrorMessage(stdlink);
          if (!!MLPutInteger32( stdlink,thestat.stepcount)) MLErrorMessage(stdlink);
          if (!!MLPutInteger32( stdlink,thestat.totitcount)) MLErrorMessage(stdlink); 
          if (!!MLPutInteger32( stdlink,thestat.maxitcount)) MLErrorMessage(stdlink); 
          if (!!MLPutInteger32( stdlink,thestat.itzero)) MLErrorMessage(stdlink); 
          if (!!MLPutInteger32( stdlink,thestat.fcn)) MLErrorMessage(stdlink); 
          if (!MLPutReal64( stdlink,(float) (clock1 - clock0)/CLOCKS_PER_SEC)) MLErrorMessage(stdlink);
          if (!MLPutReal64( stdlink,(float) (wtime1 - wtime0))) MLErrorMessage(stdlink); 
          if (!MLPutInteger32( stdlink,thestat.convergence)) MLErrorMessage(stdlink); 
          if (!MLPutInteger32( stdlink,thestat.nout)) MLErrorMessage(stdlink);     
          if(! MLEndPacket(stdlink)) MLErrorMessage(stdlink);
          if(! MLFlush(stdlink))  MLErrorMessage(stdlink); 

    break;

    case 21: case 22:

          aux[0]=thestat.totitcount;
          aux[1]=thestat2.totitcount;

          if (!MLPutFunction(stdlink, "List",12)) MLErrorMessage(stdlink);
          if (!MLPutReal64( stdlink,thestat.MaxDE)) MLErrorMessage(stdlink); 
	  if (!MLPutByteString( stdlink, (const unsigned char *)u.uu,qluz)) MLErrorMessage(stdlink);
	  if (!MLPutByteString( stdlink, (const unsigned char *)u.ee,elluz)) MLErrorMessage(stdlink);
          if (!!MLPutInteger32( stdlink,thestat.stepcount)) MLErrorMessage(stdlink);
          if (!!MLPutInteger32List (stdlink,aux,2)) MLErrorMessage(stdlink); 
          if (!!MLPutInteger32( stdlink,thestat.maxitcount)) MLErrorMessage(stdlink); 
          if (!!MLPutInteger32( stdlink,thestat.itzero)) MLErrorMessage(stdlink);
          if (!!MLPutInteger32( stdlink,thestat.fcn)) MLErrorMessage(stdlink); 
          if (!MLPutReal64( stdlink,(float) (clock1 - clock0)/CLOCKS_PER_SEC)) MLErrorMessage(stdlink); 
          if (!MLPutReal64( stdlink,(float) (wtime1 - wtime0))) MLErrorMessage(stdlink); 
          if (!MLPutInteger32( stdlink,thestat.convergence)) MLErrorMessage(stdlink);   
          if (!MLPutInteger32( stdlink,thestat.nout)) MLErrorMessage(stdlink);  
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
    free(gsmethod.b); 
    free(gsmethod.hb);      
    free(gsmethod.c); 
    free(gsmethod.hc); 
    free(gsmethod.nu); 
    free(gsmethod.orderedindices);

    free(thestat.z); 
    free(thestat.li);


}


int main(int argc, char *argv[])
{
     return MLMain(argc, argv);
}


