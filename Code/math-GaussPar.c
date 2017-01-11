/*----------------------------------------------------------------------------*/
/*									      */
/*                                math-GaussPar.c			      */
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

/* Global variables */

int thread_count;

void mathGaussPar (int neq, int n, int ns, double t0, double tend,
                   double *u0,long ulen, double *e0,long elen, 
                   val_type h, double *rpar,long rlen, int *ipar,long ilen,
                   int approximation, int threads,int algorithm1,int algorithm2,
                   int rdigits1,int rdigits2,
                   const char *myfilename1,const char *myfilename2,
                   int sampling,int codfun)

{
/*------ declarations --------------------------------------------------------*/

     int i;
     solution u; 
     solver_stat thestat;
     gauss_method gsmethod;
     ode_sys system;
     parameters params;
     toptions options;
     int aux1[2],aux2[2];
 
     u.uu = (val_type *)malloc(MAXNEQ*sizeof(val_type));
     u.ee = (val_type *)malloc(MAXNEQ*sizeof(val_type));

     options.rtol=malloc(MAXNEQ*sizeof(val_type));
     options.atol=malloc(MAXNEQ*sizeof(val_type));

     clock_t clock0, clock1; 
     time_t  wtime0,wtime1;
 
     options.rdigits=rdigits1;       
	    
/*----- Second integrations variables ----------------------------------------*/

     solution u2;
     toptions options2;
     parameters params2;
     ode_sys system2;
     solver_stat thestat2;

     params2.rpar =(val_type *)malloc(MAXPARAM*sizeof(val_type));
     params2.ipar =(int *)malloc(MAXPARAM*sizeof(int));

     u2.uu = (val_type *)malloc(MAXNEQ*sizeof(val_type));
     u2.ee = (val_type *)malloc(MAXNEQ*sizeof(val_type));

     options2.rtol=malloc(MAXNEQ*sizeof(val_type));
     options2.atol=malloc(MAXNEQ*sizeof(val_type));

     options2.rdigits=rdigits2;

/* ----------- implementation  -----------------------------------------------*/     

     params.rpar =(val_type *)malloc(MAXPARAM*sizeof(val_type));
     params.ipar =(int *)malloc(MAXPARAM*sizeof(int));
  
     thread_count=threads;
     options.t0=t0;
     options.t1=tend;
     options.h=h;
     options.algorithm=algorithm1;
     options.sampling=sampling;
     options.approximation=approximation;
     gsmethod.ns=ns;
     system.neq=neq;
     system.n=n;

/* -----Second integration initialization ------------------------------------*/

     options2.t0=options.t0;
     options2.t1=options.t1;
     options2.h = options.h;        			 
     options2.sampling=options.sampling;
     options2.approximation=options.approximation; 
     options2.algorithm=algorithm2;   		        
     options2.approximation=options.approximation;     	          

     system2.neq=system.neq;
     system2.n=system.n;
   
     strncpy(thestat2.filename, myfilename2,STRMAX); 
     InitStat(&system,&gsmethod,&thestat2);	

/* ----------- execution  ----------------------------------------------------*/

     for (i=0; i<rlen; i++) 
     {
           params.rpar[i]=rpar[i];
           params2.rpar[i]=rpar[i];
     }
 
     for (i=0; i<ilen; i++) 
     {
           params.ipar[i]=ipar[i];
           params2.ipar[i]=ipar[i];
     }

     switch (codfun)
     { 
     case 1: 
           system.f = Ode1;
           system.ham= Ham1;
           system2.f=system.f;
           system2.ham=system.ham;
     break;

     case 2: 
           system.f = Ode2;
           system.ham= Ham2;
           system2.f = system.f;
           system2.ham= system.ham; 
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
     system2.params.rpar=&params2.rpar[0];
     system2.params.ipar=&params2.ipar[0];

     for (i=0; i<system.neq; i++)        
     {
           options.rtol[i]=RTOL;
           options.atol[i]=ATOL;
           options2.rtol[i]=RTOL;
           options2.atol[i]=ATOL;
     }

     if (options.rdigits>0) options.mrdigits=pow(2,options.rdigits);
     if (options2.rdigits>0) options2.mrdigits=pow(2,options2.rdigits);
  
     GaussCoefficients(DIR_MATH,&gsmethod,&options);

     for (i=0; i<neq;i++)
     { 
           u.uu[i]=u0[i];
           u.ee[i]=e0[i];		
           u2.uu[i]=u0[i];
           u2.ee[i]=e0[i];				
     }

     strncpy(thestat.filename, myfilename1,STRMAX); 
     InitStat(&system,&gsmethod,&thestat);

     wtime0= time(NULL);
     clock0= clock();

/* ---- First integrations options -------------------------------------------*/

     switch (options.algorithm)
     {   
     case  1:  /* Jacobi*/
	   options.iteration[0]=It_Jacobi;
           system.cod[0]=0;
     break;  
     
     case  11:  /* Seidel */
           options.iteration[0]=It_Seidel;
           system.cod[0]=1;
           system.cod[1]=2;                          
     break;  

     default:
           printf("Error: incorrect algorithm\n");
     break;        
     } 

/*----- Second integrations options ------------------------------------------*/

     switch (options2.algorithm)
     {   
     case  1:  /* Jacobi*/
           options2.iteration[0]=It_Jacobi;
           system2.cod[0]=0;
     break;  
  
     case  11:  /* Seidel*/
           options2.iteration[0]=It_Seidel;
           system2.cod[0]=1;
           system2.cod[1]=2;                
     break;  

     default:
           printf("Error: incorrect algorithm\n");
     break;
     } 


#    pragma omp parallel num_threads(thread_count)
     {
#          pragma omp sections
           {
#                pragma omp section
                 {
                 /* Exec integration 1*/
                 RKG (&gsmethod,&u,&system,&options,Fixed_point_it,&thestat);                             
                 }
#                pragma omp section
                 {
                 /* Exec integration 2*/
                 RKG (&gsmethod,&u2,&system2,&options2,Fixed_point_it,&thestat2);                             
                 }
            }
      }


     clock1=clock();
     wtime1= time(NULL);

/* ----- Note !!! Return the first integration results -----------------------*/

     aux1[0]=thestat.totitcount;
     aux1[1]=thestat2.totitcount;
     aux2[0]=thestat.itzero;
     aux2[1]=thestat2.itzero;


     if (!MLPutFunction(stdlink, "List",12)) MLErrorMessage(stdlink);
     if (!MLPutReal64( stdlink,thestat.MaxDE)) MLErrorMessage(stdlink); 
     if (!MLPutReal64List( stdlink, u.uu, neq))  MLErrorMessage(stdlink); 
     if (!MLPutReal64List( stdlink, u.ee, neq))  MLErrorMessage(stdlink); 
     if (!!MLPutInteger32( stdlink,thestat.stepcount)) MLErrorMessage(stdlink);
     if (!!MLPutInteger32List( stdlink,aux1,2)) MLErrorMessage(stdlink); 
     if (!!MLPutInteger32( stdlink,thestat.maxitcount)) MLErrorMessage(stdlink); 
     if (!!MLPutInteger32List( stdlink,aux2,2)) MLErrorMessage(stdlink); 
     if (!!MLPutInteger32( stdlink,thestat.fcn)) MLErrorMessage(stdlink); 
     if (!MLPutReal64( stdlink,(float) (clock1 - clock0)/CLOCKS_PER_SEC)) 
         MLErrorMessage(stdlink); 
     if (!MLPutReal64( stdlink,(float) (wtime1 - wtime0)))
         MLErrorMessage(stdlink); 
     if (!MLPutInteger32( stdlink,thestat.convergence)) MLErrorMessage(stdlink);  
     if (!MLPutInteger32( stdlink,thestat.nout)) MLErrorMessage(stdlink);    
     if(! MLEndPacket(stdlink)) MLErrorMessage(stdlink);
     if(! MLFlush(stdlink))  MLErrorMessage(stdlink); 

     free(params.rpar); free(params.ipar);
     free(params2.rpar); free(params2.ipar);
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


