/*----------------------------------------------------------------------------*/
/*									      */
/*    GaussCommon.c							      */
/*									      */
/*	Functions: 							      */
/*       print_u():							      */
/*	 InitStat():							      */
/*	 NormalizedDistance():						      */
/*       StatYinit():							      */
/*	 RemoveDigits():						      */
/*	 Default_Stage_init();						      */
/*	 Interpolated_Stage_init();					      */
/*	 StopCriterion();						      */
/*	 Fixed_point_Step():						      */
/*	 General_FP_It():						      */
/*	 Partitioned_FP_It():						      */
/*	 MyOutput():							      */
/*	 CompensatedSummation():					      */
/*	 IRKFP():							      */
/*       								      */
/* ---------------------------------------------------------------------------*/

#include <Common-IRKFP.h>
#include <float.h> 

/******************************************************************************/
/* 					   				      */
/* print_u 		         	  				      */
/*                                        				      */
/*									      */
/******************************************************************************/

void print_u (const int neq, const val_type *u)
{
     int i;

     for (i = 0;i<neq;i++) 
     {
           if (i<neq-1)
               printf("%.20lg,", u[i]);
           else
           printf("%.20lg\n", u[i]);
     }     

}


/******************************************************************************/
/* 					   				      */
/* InitStat 		         	  				      */
/*                                        				      */
/*									      */
/******************************************************************************/
      
void InitStat (const ode_sys *system,const gauss_method *gsmethod, 
               solver_stat *thestatptr)
{

     int i,is,neq,ns;

     ns=gsmethod->ns;
     neq=system->neq;

     thestatptr->laststep = false;    
     thestatptr->stepcount = 0;
     thestatptr->itcount=0;
     thestatptr->totitcount=0;
     thestatptr->totitcountzero=0;
     thestatptr->maxitcount=0;
     thestatptr->itzero=0;
     thestatptr->fcn=0;
     thestatptr->convergence=SUCCESS;
     thestatptr->nout=0;
     thestatptr->MaxDE=0.;
      
     thestatptr->z = (val_type *)malloc(neq*(ns)*sizeof(val_type));
     thestatptr->li = (val_type *)malloc(neq*ns*sizeof(val_type));
     thestatptr->fz = (val_type *)malloc(neq*ns*sizeof(val_type));
     thestatptr->zit0 = (val_type *)malloc(neq*ns*sizeof(val_type));
     thestatptr->initqlty = (int *)malloc(neq*ns*sizeof(int));

     for (is=0; is<ns; is++)
          for (i=0; i<neq; i++)
          {
               thestatptr->li[is*neq+i]=0.;
               thestatptr->fz[is*neq+i]=0.;
               thestatptr->initqlty[is*neq+i]=0;
          }

     return;

}


/******************************************************************************/
/* 					   				      */
/* NormalizedDistance: 	         	  				      */
/*                                        				      */
/*									      */
/******************************************************************************/

val_type NormalizedDistance ( const int neq, const int ns,
                              const toptions *options, const val_type *z,
                              const val_type *zold)
{


/*---------------- Declarations ----------------------------------------------*/

     int i,is,ix;
     val_type maxi,mayi,relerrors;
     val_type maxz,maxzold;

/* --------------- Implementation --------------------------------------------*/
   	
     maxi=0.;

     for (i=0; i<neq; i++)
     {  
          maxz=0.;
          maxzold=0.;
            
          for (is=0; is<ns;is++)
          {
               ix=neq*is+i;
               maxz=FMAX(maxz,FABS(z[ix]));
               maxzold=FMAX(maxzold,FABS(zold[ix]));
          }

          mayi=(maxz+maxzold)/2;
          relerrors=0.;
      
          for (is=0; is<ns; is++)
          {
               ix=neq*is+i;
               relerrors=FMAX(relerrors,
                              FABS(z[ix]-zold[ix])/
                                   (mayi*options->rtol[i]+options->atol[i]));
          }

          maxi=FMAX(maxi,relerrors);

     }

     return maxi;

}


/******************************************************************************/
/* 					   				      */
/* StatYinit: 		         	  				      */
/*                                        				      */
/*									      */
/******************************************************************************/

int StatYinit (const ode_sys *system,const gauss_method *method,
               solver_stat *thestatptr)
{


/* ---------- First initializations ------------------------------------------*/

     int neq,ns;

     neq=system->neq;
     ns=method->ns;

/*------------- declarations -------------------------------------------------*/ 

     int i,is,in;
     val_type *z,*zit0;
     int *initqlty;
     int prec;

/* ----------- implementation  -----------------------------------------------*/

     z=thestatptr->z;
     zit0=thestatptr->zit0;
     initqlty=thestatptr->initqlty;
     prec=16;

     for (is = 0; is<ns; is++)
           for (i = 0; i<neq; i++)
           {                      
                 in=neq*is+i;
                 if (z[in]==zit0[in]) {initqlty[in]+=prec;}
                       else 
                       {
                        initqlty[in]+=(int)round(fabs(log10(fabs(z[in]-zit0[in])/fabs(zit0[in]))));
                       }
           }
                    
     return(0);

}



/******************************************************************************/
/* 					   				      */
/* RemoveDigits		         	  				      */
/*                                        				      */
/*									      */
/******************************************************************************/

void RemoveDigits (val_type *x,const int m)
{

/* ---------- First initializations ------------------------------------------*/

/*---------------- Declarations ----------------------------------------------*/

     val_type aux,mx,mxx;

/* --------------- Implementation --------------------------------------------*/

     aux=*x;
     mx=m*aux;
     mxx=mx+aux;
     *x=mxx-mx;

     return;

}

/******************************************************************************/
/* 					   				      */
/* Default_Stage_init:	         	  				      */
/*                                        				      */
/*									      */
/******************************************************************************/

void Default_Stage_init (const solution *u,  val_type *z, const ode_sys *system,
             const gauss_method *method,solver_stat *thestatptr,
             const toptions *options)
{


/* ---------- First initializations ------------------------------------------*/

     int neq,ns;

     neq=system->neq;
     ns=method->ns;

/*------------- declarations -------------------------------------------------*/ 

     int i,is;
     val_type *zit0;

/* ----------- implementation  -----------------------------------------------*/
       

          zit0=thestatptr->zit0;
          for (is = 0; is<ns; is++)
               for (i = 0; i<neq; i++)
               {
                    z[is*neq+i]=u->uu[i];
                    zit0[is*neq+i]=u->uu[i];   
               }
    
 
     return;

}

/******************************************************************************/
/* 					   				      */
/* Interpolated_Stage_init:         	  				      */
/*                                        				      */
/*									      */
/******************************************************************************/

void Interpolated_Stage_init (const solution *u,  val_type *z, const ode_sys *system,
             const gauss_method *method,solver_stat *thestatptr,
             const toptions *options)
{


/* ---------- First initializations ------------------------------------------*/

     int neq,ns;

     neq=system->neq;
     ns=method->ns;

/*------------- declarations -------------------------------------------------*/ 

     int i,is,in,js,jsn;
     val_type *coef,*li,*zit0;
     val_type sum;

/* ----------- implementation  -----------------------------------------------*/
        
   
          li=thestatptr->li;
          zit0=thestatptr->zit0;
          coef=method->nu;     
  
          for (is = 0; is<ns; is++)
          {
                for (i = 0; i<neq; i++)
                {  
                     in=neq*is+i;
                     sum=0.;

                     for (js = 0; js<ns; js++)
                     {
                          jsn=ns*is+js;
                          sum+=li[neq*js+i]*(coef[jsn]);
                     }

                     z[in]=u->uu[i]+sum;
                     zit0[in]=z[in];
                     
                }
          }  
 
     return;

}


/******************************************************************************/
/* 									      */
/*   StopCriterion  						              */
/* 									      */
/* 									      */
/******************************************************************************/

void StopCriterion    (const ode_sys *system, const gauss_method *method,
                       int *D0,bool *cont,val_type *DMin, 
                       const val_type *Y, const val_type *Yold)
{
/*----------------  First: initialization ------------------------------------*/
     int neq,ns;

     neq=system->neq;
     ns=method->ns;

/*------ declarations --------------------------------------------------------*/
     int i,is,isn;
     val_type dY;

/* ----------- implementation  -----------------------------------------------*/ 
 

     bool plusIt;

     if (*D0<0) plusIt=false;
         else plusIt=true;


     *D0=0;
     *cont=false;


     for (is=0; is<ns; is++)
          for (i=0; i<neq; i++) 
          { 
               isn=is*neq+i;
               dY=FABS(Y[isn]-Yold[isn]);

               if (dY>0.)
               {                     
                   if (dY<DMin[isn])
                    {
                         DMin[isn]=dY;
                         *cont=true;
                    }

               }
               else
               {
                    *D0=*D0+1;
               } 
          }
 
    if (*cont==false && *D0<(ns*neq) && plusIt)
    {
          *D0=-1;
          *cont=true;
    }

     return;

}




/******************************************************************************/
/* 									      */
/*   Fixed_point_Step  							      */
/* 									      */
/* 									      */
/******************************************************************************/


void Fixed_point_Step ( const ode_sys *system, const solution *u, 
                      const val_type tn, const val_type h, 
                      const toptions *options, const gauss_method *method,
                      solver_stat *thestatptr)

{ 

/*----------------  First: initialization ------------------------------------*/
     int neq,ns;

     neq=system->neq;
     ns=method->ns;

/*------ declarations --------------------------------------------------------*/

     int i,is;
     bool iter0;
     int D0;
     val_type difftest,DMin[neq*ns];
     val_type *z;
     val_type zold[neq*ns];


/* ----------- implementation  -----------------------------------------------*/ 

     z=thestatptr->z;   

     iter0=true;
     D0=0;
     for (is=0; is<ns; is++)
          for (i=0; i<neq; i++) DMin[is*neq+i]=INF; 

     thestatptr->itcount=0;

     while (iter0 && thestatptr->itcount<MAXIT)
     { 

           for (is=0; is<ns;is++)
                for (i=0; i<neq; i++) zold[neq*is+i]=thestatptr->z[neq*is+i];

           options->iteration(system,u,tn,h,method,thestatptr); 
           StopCriterion(system,method,&D0,&iter0,DMin,z,zold);
           thestatptr->itcount++;
     }

     if (thestatptr->itcount==MAXIT)
     { 
           printf("Break: step(MAXIT)=%i\n",thestatptr->itcount);
           thestatptr->convergence=FAIL; 
     }
 
     if (D0<(ns*neq))
     {
           difftest=NormalizedDistance(neq,ns,options,z,zold);
           if (difftest>1.)
           {
                thestatptr->convergence=FAIL;
                printf("Lack of convegence of Fixed point iteration:\
                        step=%i,iteration=%i,",
                        thestatptr->stepcount,thestatptr->itcount);  
                printf("difftest=%lg\n",difftest);                   
           }
     }
     else
     {
           (thestatptr->totitcountzero)+=(thestatptr->itcount);
           thestatptr->itzero++;
     }

     return;

}


/******************************************************************************/
/*									      */
/*      General_FP_It: General Fixed Point Iteration			      */
/*									      */
/******************************************************************************/


int General_FP_It (const ode_sys *system, const solution *u, const val_type tn,
                   const val_type h, const gauss_method *method,
                   solver_stat *thestatptr)
{

#ifdef PARALLEL
     int extern thread_count;
#endif

/* ---------- First initializations ------------------------------------------*/

     int neq,ns;
     parameters *params;
     val_type *z,*li,*fz;
    
     neq=system->neq;
     ns=method->ns;
     params=&system->params;

     z=thestatptr->z;
     li=thestatptr->li;
     fz=thestatptr->fz; 

/*------ declarations --------------------------------------------------------*/

     int i,is,js,isn,in,jsn;
     val_type sum;
  
/* ----------- implementation  -----------------------------------------------*/ 

#ifdef PARALLEL
#      pragma omp parallel for num_threads(thread_count) private(isn)
#endif  
     for (is = 0; is<ns; is++)
     {
           isn=neq*is;
           params->eval=system->cod[0];
           thestatptr->fcn++;
           system->f(neq,tn+method->hc[is],&z[isn],&fz[isn],params); 
           for (i=0; i<neq; i++) li[isn+i]=fz[isn+i]*method->hb[is];		
     }

     for (is = 0; is<ns; is++)
     { 
           for (i = 0; i<neq; i++)    
           {  
                 in=neq*is+i; 
                 sum=method->m[ns*is]*li[i]+u->ee[i];                       
                 for (js =1; js<ns; js++)
                 {
                       jsn=ns*is+js;
                       sum+=method->m[jsn]*li[neq*js+i];
                 }
                 z[in]=u->uu[i]+sum;        
           }
     } 
   
     return(0);

}	

/******************************************************************************/
/*									      */
/*      Partitioned_FP_It: Partitioned Fixed Point Iteration		      */
/*									      */
/******************************************************************************/

int Partitioned_FP_It
(const ode_sys *system, const solution *u, const val_type tn,
 const val_type h, const gauss_method *method,
 solver_stat *thestatptr)

{

#ifdef PARALLEL
     int extern thread_count;
#endif

/* ---------- First initializations ------------------------------------------*/

     int neq,ns;
     parameters *params;
     val_type *z,*li,*fz;
    
     neq=system->neq;
     ns=method->ns;
     params=&system->params;

     z=thestatptr->z;
     li=thestatptr->li;
     fz=thestatptr->fz;

/*------ declarations --------------------------------------------------------*/

     int i,is,js,isn,in,jsn;
     val_type sum;
  
/* ----------- implementation  -----------------------------------------------*/ 


/* ----------- first part  ---------------------------------------------------*/ 
#ifdef PARALLEL
#      pragma omp parallel for num_threads(thread_count) private(isn)
#endif  
     for (is = 0; is<ns; is++)
     {
           isn=neq*is;
           params->eval=system->cod[0];
           thestatptr->fcn++;
           system->f(neq,tn+method->hc[is],&z[isn],&fz[isn],params);
           for (i=0; i<neq/2; i++) li[isn+i]=fz[isn+i]*method->hb[is];	
     }

     for (is = 0; is<ns; is++)
     { 
           for (i = 0; i<neq/2; i++)    
           {  
                 in=neq*is+i; 
                 sum=method->m[ns*is]*li[i]+u->ee[i];
                         
                 for (js =1; js<ns; js++)
                 {
                       jsn=ns*is+js;
                       sum+=method->m[jsn]*li[neq*js+i];
                 }

                 z[in]=u->uu[i]+sum;        
           }
     } 


/* ----------- second part  --------------------------------------------------*/ 

#ifdef PARALLEL
#      pragma omp parallel for num_threads(thread_count) private(isn)
#endif  
     for (is = 0; is<ns; is++)
     {
           isn=neq*is;
           params->eval=system->cod[1];
           system->f(neq,tn+method->hc[is],&z[isn],&fz[isn],params);
           for (i=neq/2; i<neq; i++) li[isn+i]=fz[isn+i]*method->hb[is];		
      }

     for (is = 0; is<ns; is++)
     { 
           for (i=neq/2; i<neq; i++)    
           {  
                 in=neq*is+i; 
                 sum=method->m[ns*is]*li[i]+u->ee[i];
                         
                 for (js =1; js<ns; js++)
                 {
                       jsn=ns*is+js;
                       sum+=method->m[jsn]*li[neq*js+i];
                 }

                 z[in]=u->uu[i]+sum;        
            }
     } 
    
     return(0);

}	



/******************************************************************************/
/*									      */
/*      MyOutput 							      */
/*									      */
/******************************************************************************/
 
void MyOutput(const ode_sys *system,const gauss_method *method,
               const val_type t, val_type h, const solution *u, 
               solver_stat *thestatptr,
               const parameters *params,const toptions *options,FILE *myfile)
{

/* ---------- First initializations ------------------------------------------*/

     int neq;
     neq=system->neq;


/*------ declarations --------------------------------------------------------*/
     val_type DH;
     int i,node;
  
     struct rec
     {
     val_type t;
     val_type uu[neq];
     val_type ee[neq];
     };
  
     struct rec my_record;

/* ----------- implementation  -----------------------------------------------*/ 

     node=neq;
     DH=(system->ham(node,u,params)-thestatptr->E0)/thestatptr->E0;
     if (FABS(DH)>thestatptr->MaxDE) thestatptr->MaxDE=FABS(DH);

     if (((thestatptr->stepcount % options->sampling) == 0) || (thestatptr->laststep))
     {
           thestatptr->nout++;
           printf("%i,%i, %lg\n",thestatptr->stepcount,thestatptr->itcount,DH);

           my_record.t=t;
           for (i=0; i<neq;i++) 
           {
                 my_record.uu[i]=u->uu[i];
                 my_record.ee[i]=u->ee[i];
           }
           fwrite(&my_record, sizeof(struct rec), 1, myfile);
  
     } 

     if (thestatptr->stepcount>1) StatYinit(system,method,thestatptr);

  return;

}


/******************************************************************************/
/* 									      */
/*   Compensated Summation:						      */
/* 									      */
/* 									      */
/******************************************************************************/

void CompensatedSummation (const gauss_method *gsmethod,
                           solution *u,const ode_sys *system,
                           const toptions *options,const solver_stat *thestatptr)

{
/* ---------- First initializations ------------------------------------------*/

     int neq,ns;
     
     neq=system->neq;
     ns=gsmethod->ns;

/*------ declarations --------------------------------------------------------*/

     int i,ix,is,isn;
     val_type aux,*li;
     li=thestatptr->li;
  
     val_type s0,s1,delta,eli,ee;
     val_type *fz;
     fz=thestatptr->fz;


/* ----------- implementation  -----------------------------------------------*/ 


     /* --------------- update (e_n= en+ Sum ELj) ----------------------------*/


          for (i = 0; i<neq; i++)
          {
              eli= fz[i]*gsmethod->hb[0]-li[i]; 

              for (is=1; is<ns; is++)
              {
                      isn=neq*is+i;
                      aux=fz[isn]*gsmethod->hb[is]-li[isn];
                      eli+=aux;                      
              }

              u->ee[i]+=eli;
       
          }


     /* ----------------yn+1=yn+(Sum Li+e_n) 21-07-2016-----------------------*/

          for (i = 0; i<neq; i++)
          { 

               s0=u->uu[i];
               ee=u->ee[i];
		
               for (ix =0 ; ix<ns; ix++)
               {
                    is=gsmethod->orderedindices[ix];
                    isn=neq*is+i;                
                    delta=li[isn]+ee;
                    if (options->rdigits>0) 
                        RemoveDigits(&delta,options->mrdigits);
                    s1=s0;
                    s0=s1+delta;
                    aux=s1-s0;
                    ee=aux+delta;
               } 

               u->uu[i]=s0;
               u->ee[i]=ee;
      	       	       
          }


     return;
}




/******************************************************************************/
/* 									      */
/*   IRKFP  								      */
/* 									      */
/* 									      */
/******************************************************************************/

void IRKFP
(val_type t0, val_type t1, val_type h,
 const gauss_method *gsmethod, solution *u,
 const ode_sys *system, toptions *options,
 solver_stat *thestatptr)

{

/* ---------- First initializations ------------------------------------------*/

     int neq;
     parameters *params;
     
     neq=system->neq;
     params=&system->params;

/*------ declarations --------------------------------------------------------*/

     FILE *myfile;
  
     int istep,nstep;
     val_type tn;             
     val_type *z;

     z=thestatptr->z;

/* ----------- implementation  -----------------------------------------------*/ 

     if (options->TheOutput != 0)
           myfile = fopen(options->filename,"wb"); 

/* ----------- initial energy (E0)   -----------------------------------------*/

     thestatptr->E0=system->ham(neq,u,params);
     printf("Initial energy=%lg\n", thestatptr->E0);

     tn=t0;
     nstep=(t1-t0)/h; 

     if (options->TheOutput != 0)
          options->TheOutput(system,gsmethod,tn,h,u,thestatptr,params,options,myfile);

     for(istep=0; istep<nstep; istep++) 
     {     
          options->StageInitFn (u,z,system,gsmethod,thestatptr,options);      
          Fixed_point_Step (system,u,tn,h,options,gsmethod,thestatptr);
          if (thestatptr->convergence==FAIL)
          { 
               printf("Stop Fail. step=%i\n", istep);
               nstep=istep;
          } 	      
          else
          {
               CompensatedSummation (gsmethod,u,system,options,thestatptr);

               tn=(istep+1)*h;

               thestatptr->stepcount++;		
               (thestatptr->totitcount)+=(thestatptr->itcount); 
               if ((thestatptr->itcount)>(thestatptr->maxitcount))
                  (thestatptr->maxitcount)=(thestatptr->itcount);

               if (options->TheOutput != 0)
                   options->TheOutput(system,gsmethod,tn,h,u,thestatptr,params,options,myfile);

               if (tn+h>=t1)
               { 
                    h=t1-tn;
                    thestatptr->laststep=true;
               } 
          }
     }

     if (options->TheOutput != 0) fclose(myfile);
     return;

}








