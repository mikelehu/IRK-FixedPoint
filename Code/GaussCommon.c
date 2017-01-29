/*----------------------------------------------------------------------------*/
/*									      */
/*    GaussCommon.c							      */
/*									      */
/*	Functions: 							      */
/*	 InitStat():							      */
/*	 NormalizedDistance():						      */
/*       statYinit():							      */
/*	 RemoveDigitsFcn():						      */
/*	 Yi_init();							      */
/*	 StopCriterion();						      */
/*	 Fixed_point_it():						      */
/*	 It_Jacobi():							      */
/*	 TheOutput():							      */
/*	 TheOutput2():							      */
/*	 RKG():								      */
/*	 RKG2():							      */
/*	 CompensatedSummation():					      */
/*       select_gauss();						      */
/*       select_odefun();						      */
/* ---------------------------------------------------------------------------*/

#include <GaussCommon.h>
#include <quadmath.h>
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
/* statzinit: 		         	  				      */
/*                                        				      */
/*									      */
/******************************************************************************/

int statYinit (const ode_sys *system,const gauss_method *method,
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
/* RemoveDigitsFcn	         	  				      */
/*                                        				      */
/*									      */
/******************************************************************************/

void RemoveDigitsFcn (val_type *x,const int m)
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
/* Yi_init: 		         	  				      */
/*                                        				      */
/*									      */
/******************************************************************************/

int Yi_init (const solution *u,  val_type *z, const ode_sys *system,
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

 
     switch (options->approximation)
     {

     case 0:           

          zit0=thestatptr->zit0;
          for (is = 0; is<ns; is++)
               for (i = 0; i<neq; i++)
               {
                    z[is*neq+i]=u->uu[i];
                    zit0[is*neq+i]=u->uu[i];   
               }
     break;

     case 1:          
   
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
     break;
          
     default:
          printf("approximation=%i error\n", options->approximation);
     break;

     }
 
     return(0);

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
/*   Fixed_point_it  							      */
/* 									      */
/* 									      */
/******************************************************************************/


void Fixed_point_it ( const ode_sys *system, const solution *u, 
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
  
     for (is=0; is<ns; is++)
          for (i=0; i<neq; i++) DMin[is*neq+i]=INF; 

     thestatptr->itcount=0;

     while (iter0 && thestatptr->itcount<MAXIT)
     { 

           for (is=0; is<ns;is++)
                for (i=0; i<neq; i++) zold[neq*is+i]=thestatptr->z[neq*is+i];

           options->iteration[0](system,u,tn,h,method,thestatptr); 
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
/*      It_Jacobi							      */
/*									      */
/******************************************************************************/


int It_Jacobi (const ode_sys *system, const solution *u, const val_type tn,
               const val_type h, const gauss_method *method,
               solver_stat *thestatptr)
{

     int extern thread_count;

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
/*      TheOutput (RKG)							      */
/*									      */
/******************************************************************************/
 
void TheOutput(const ode_sys *system,const gauss_method *method,
               const val_type t, const solution *u,
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

     if (thestatptr->stepcount>1) statYinit(system,method,thestatptr);

  return;

}

/******************************************************************************/
/*									      */
/*      TheOutput2 (RKG2)						      */
/*									      */
/******************************************************************************/

void TheOutput2(const ode_sys *system, const gauss_method *method,
                const val_type t, const solution *u,
                const solution *u2,solver_stat *thestatptr,
                const parameters *params, const toptions *options,FILE *myfile)
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
     val_type est[neq];
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
           my_record.est[i]=(u->uu[i]-u2->uu[i])+(u->ee[i]-u2->ee[i]);}
           fwrite(&my_record, sizeof(struct rec), 1, myfile);  

     } 


     if (thestatptr->stepcount>1) statYinit(system,method,thestatptr);

     return;

}


/******************************************************************************/
/* 									      */
/*   Compensated Summation:						      */
/* 									      */
/* 									      */
/******************************************************************************/

void CompensatedSummation (const gauss_method *gsmethod,
                           val_type *u0,solution *u,const ode_sys *system,
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

               u0[i]=u->uu[i];
               s0=u->uu[i];
               ee=u->ee[i];
		
               for (ix =0 ; ix<ns; ix++)
               {
                    is=gsmethod->orderedindices[ix];
                    isn=neq*is+i;                
                    delta=li[isn]+ee;
                    if (options->rdigits>0) 
                        RemoveDigitsFcn(&delta,options->mrdigits);
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
/*   RGK:  								      */
/* 									      */
/* 									      */
/******************************************************************************/

void RKG
(const gauss_method *gsmethod, solution *u,
 const ode_sys *system, toptions *options,
 void RKG_Step (), solver_stat *thestatptr)

{

/* ---------- First initializations ------------------------------------------*/

     int neq;
     parameters *params;
     
     neq=system->neq;
     params=&system->params;

/*------ declarations --------------------------------------------------------*/

     FILE *myfile;
  
     int istep,nstep;
     val_type u0[neq];
     val_type tn;             
     val_type *z;

     z=thestatptr->z;

/* ----------- implementation  -----------------------------------------------*/ 

     myfile = fopen(thestatptr->filename,"wb"); 

/* ----------- initial energi (E0)   -----------------------------------------*/

     thestatptr->E0=system->ham(neq,u,params);
     printf("Initial energy=%lg\n", thestatptr->E0);

     tn=options->t0;
     nstep=((options->t1)-(options->t0))/(options->h); 

#ifdef IOUT
     TheOutput(system,gsmethod,tn,u,thestatptr,params,options,myfile);
#endif 

     for(istep=0; istep<nstep; istep++) 
     {     
          Yi_init (u,z,system,gsmethod,thestatptr,options);      
          RKG_Step (system,u,tn,options->h,options,gsmethod,thestatptr);
          if (thestatptr->convergence==FAIL)
          { 
               printf("Stop Fail. step=%i\n", istep);
               nstep=istep;
          } 	      
          else
          {
               CompensatedSummation (gsmethod,u0,u,system,options,thestatptr);

               tn=(istep+1)*(options->h);

               thestatptr->stepcount++;		
               (thestatptr->totitcount)+=(thestatptr->itcount); 
               if ((thestatptr->itcount)>(thestatptr->maxitcount))
                  (thestatptr->maxitcount)=(thestatptr->itcount);

#ifdef IOUT
               TheOutput(system,gsmethod,tn,u,thestatptr,params,options,myfile);
#endif 


               if ((tn+options->h)>=options->t1)
               { 
                    options->h=options->t1-tn;
                    thestatptr->laststep=true;
               } 
          }
     }


     fclose(myfile);
     return;

}


/******************************************************************************/
/* 								       	      */
/*   RKG2:  double integration (sequentially).  			      */
/*        								      */
/* 									      */
/******************************************************************************/

void RKG2 
(const gauss_method *gsmethod, const gauss_method *gsmethod2,
 solution *u,solution *u2,
 const ode_sys *system, toptions *options, toptions *options2,
 void RKG_Step (), solver_stat *thestatptr,solver_stat *thestatptr2)
 
{

/* ---------- First initializations ------------------------------------------*/

     int neq,ns;
     parameters *params;

     neq=system->neq;
     ns=gsmethod->ns;
     params=&system->params;

/*------ declarations --------------------------------------------------------*/

     FILE *myfile;
  
     int i,is,ksw;
     int istep,nstep;
     val_type tn;             

     val_type u0[neq];  
     val_type *z,*z2;
     bool initwithfirstintegration;         

     z=thestatptr->z;
     z2=thestatptr2->z;


/* ----------- implementation  -----------------------------------------------*/ 

     myfile = fopen(thestatptr->filename,"wb"); 

/* ----------- initial energi (E0)   -----------------------------------------*/

     thestatptr->E0=system->ham(neq,u,params);
     printf("Initial energy=%lg\n", thestatptr->E0);
 
     tn=options->t0;
     nstep=((options->t1)-(options->t0))/(options->h);
     if (nstep*(options->h)+(options->t0)<options->t1) 
           nstep++;    
  
     initwithfirstintegration=true;
     ksw=0;


#ifdef IOUT
     TheOutput2(system,gsmethod,tn,u,u2,thestatptr,params,options,myfile);
#endif 

     for(istep=0; istep<nstep; istep++) 
     {     

          /* ---------------- First integration ------------------------------*/

          Yi_init (u,z,system,gsmethod,thestatptr,options);   
          RKG_Step (system,u,tn,options->h,options,gsmethod,thestatptr);
 
          switch (thestatptr->convergence)  
          {
          case FAIL:
               printf("Stop Fail. step=%i\n", istep);
               nstep=istep;
          break;

          default: 	
               CompensatedSummation (gsmethod,u0,u,system,options,thestatptr);

          /* ---------------- Second integration -----------------------------*/

               if (initwithfirstintegration==false)
               {
                     Yi_init (u2,z2,system,gsmethod2,thestatptr2,options2);
               }  
               else
               {    for (is=0; is<ns; is++)
                         for (i=0; i<neq; i++)
                             {z2[is*neq+i]=z[is*neq+i]+(u2->uu[i]-u0[i]);}
               }

               RKG_Step (system,u2,tn,options2->h,options2,gsmethod2,thestatptr2);
          
               switch (thestatptr2->convergence)
               { 
               case FAIL:
                    printf("Stop Fail. step=%i\n", istep);
                    nstep=istep;
               break;
    
               default:
                    CompensatedSummation (gsmethod2,u0,u2,system,options2,thestatptr2);

                    if (NormalizedDistance(neq,1,options,u->uu,u2->uu)>ESTERRTHRS)     
                    { 
                          printf("Stop integration. Estimated error is high=%i\n", istep);
                          nstep=istep;
                    };
     
                    if (initwithfirstintegration==true)
                    {
                           if ((thestatptr2->itcount) > thestatptr->itcount)
                           {
                                 ksw++;
                                 if (ksw>MAXKSW) {initwithfirstintegration=false;}
                           }
                           else
                           {
                                 ksw=0;
                           }
                    }

                    tn=(istep+1)*(options->h);

                    thestatptr->stepcount++;		
                    (thestatptr->totitcount)+=(thestatptr->itcount); 
                    if ((thestatptr->itcount)>(thestatptr->maxitcount))
                       {(thestatptr->maxitcount)=(thestatptr->itcount);}

                    thestatptr2->stepcount++;		
                    (thestatptr2->totitcount)+=(thestatptr2->itcount); 
                    if ((thestatptr2->itcount)>(thestatptr2->maxitcount))
                       {(thestatptr2->maxitcount)=(thestatptr2->itcount);}

#ifdef IOUT
                    TheOutput2(system,gsmethod,tn,u,u2,thestatptr,params,options,myfile);
#endif 

                    if ((tn+options->h)>=options->t1)
                    { 
                          options->h=options->t1-tn;
                          thestatptr->laststep=true;
                    } 

               break;
               } 

          break;
          }
      }

      fclose(myfile);

      return;

}



/******************************************************************************/
/* 									      */
/*   select_gauss				 			      */
/*        								      */
/* 									      */
/******************************************************************************/

void select_gauss
(gauss_method *gsmethodptr, gauss_method *gsmethod2ptr,
 solution *uptr,solution *u2ptr,
 ode_sys *systemptr,toptions *optionsptr,toptions *options2ptr,
 solver_stat *thestatptr,solver_stat *thestat2ptr)

{

/* ---------- First initializations ------------------------------------------*/


/*------ declarations --------------------------------------------------------*/

    int i;


/* ----------- implementation  -----------------------------------------------*/

    switch (optionsptr->algorithm)
    {  

   /* RKG: gure inplementazioa */

    case  1:  /* Standard fixed point iteration */
          optionsptr->iteration[0]=It_Jacobi;
          systemptr->cod[0]=0;
          RKG (gsmethodptr,uptr,systemptr,optionsptr,Fixed_point_it,thestatptr);                             
    break;  

    case  2: /* RKG2: sequential execution */
          options2ptr->t0=optionsptr->t0;
  	  options2ptr->t1=optionsptr->t1;
 	  options2ptr->h=optionsptr->h;
 	  options2ptr->algorithm=optionsptr->algorithm;
 	  options2ptr->sampling=optionsptr->sampling;
 	  options2ptr->approximation=optionsptr->approximation;

          for (i=0; i<systemptr->neq; i++)
          {
               options2ptr->rtol[i]=optionsptr->rtol[i];
               options2ptr->atol[i]=optionsptr->atol[i];
          }

          
          for (i=0; i<systemptr->neq;i++)
   	  {
                u2ptr->uu[i]=uptr->uu[i];
    		u2ptr->ee[i]=uptr->ee[i];					
          }

          optionsptr->iteration[0]=It_Jacobi;
          options2ptr->iteration[0]=It_Jacobi;
          systemptr->cod[0]=0;
          gsmethod2ptr->ns=gsmethodptr->ns;

 	  InitStat(systemptr,gsmethod2ptr,thestat2ptr);
          RKG2 (gsmethodptr,gsmethod2ptr,uptr,u2ptr,systemptr,optionsptr,options2ptr,
                Fixed_point_it,thestatptr,thestat2ptr);  

    break;

                                  
    default:
          printf("Error: incorrect algorithm\n");
    break;

         
    } 
}

/******************************************************************************/
/* 									      */
/*   select_odefun				 			      */
/*        								      */
/* 									      */
/******************************************************************************/

void select_odefun (const int codfun, ode_sys *system)
{
     switch (codfun)
     { 
     case 1: 
           system->f = Ode1;
           system->ham= Ham1;                   
     break;

     case 2: 
           system->f = Ode2;
           system->ham= Ham2;
     break;

     case 3: 
           system->f = Ode3;
           system->ham= Ham3;               
     break;

     case 4: 
           system->f = Ode4;
           system->ham= Ham4;                  
     break;

     case 5: 
           system->f = Ode5;
           system->ham= Ham5;               
     break;

     case 6: 
           system->f = Ode6;
           system->ham= Ham6;               
     break;

     case 7: 
           system->f = Ode7;
           system->ham= Ham7;                   
     break;

     case 8: 
           system->f = Ode8;
           system->ham= Ham8;                  
     break;

     case 9: 
           system->f = Ode9;
           system->ham= Ham9;                 
     break;

     case 10: 
           system->f = Ode10;
           system->ham= Ham10;                   
     break;


     default:
           printf("error. codfun\n");
     break;
     
     }

     return;
}


