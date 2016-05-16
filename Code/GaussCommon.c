/*------------------------------------------------------------------------------*/
/*										*/
/*    GaussCommon.c								*/
/*										*/
/*	Functions: 								*/
/*	 print_u():								*/
/*	 InitStat():								*/
/*	 NormalizedDistance():							*/
/*	 RemoveDigitsFcn():							*/
/*	 Yi_init();								*/
/*	 UpdateDMin();								*/
/*	 Fixed_point_it():							*/
/*	 It_Jacobi():								*/
/*	 It_Seidel():								*/
/*	 TheOutput():								*/
/*	 TheOutput2():								*/
/*	 RKG():									*/
/*	 RKG2():								*/
/*       select_gauss();							*/
/* -----------------------------------------------------------------------------*/

#include <GaussCommon.h>
#include <quadmath.h>

/************************************************************************************/
/* 					   					    */
/* print_u 		         	  					    */
/*                                        					    */
/*										    */
/************************************************************************************/


void print_u (int neq,val_type *u)
{
     int i;

#    if PREC ==2  //QUADRUPLEPRECISION
     int n;
     int width = 46;
     char buf[128];
#    endif

#if PREC ==2  //QUADRUPLEPRECISION
           printf("u:");
           for (i=0; i<neq; i++)
           { 
               n = quadmath_snprintf(buf, sizeof buf, "%+-#*.30Qe", width, u[i]);
               if (i<neq-1)
               {  if ((size_t) n < sizeof buf) printf("%s,",buf);}
               else
               {  if ((size_t) n < sizeof buf) printf("%s\n,",buf);} 
           }        

#else  // DOUBLEPRECISION, FLOAT    
           for (i = 0;i<neq;i++) 
           {
               if (i<neq-1)
                  printf("%.20lg,", u[i]);
               else
                  printf("%.20lg\n", u[i]);
           }     
#endif 


}


/************************************************************************************/
/* 					   					    */
/* InitStat 		         	  					    */
/*                                        					    */
/*										    */
/************************************************************************************/
      
void InitStat (ode_sys *system,gauss_method *gsmethod, solver_stat *thestatptr)
{

     int i,is,neq,ns;

     ns=gsmethod->ns;
     neq=system->neq;

     thestatptr->laststep = false;    
     thestatptr->stepcount = 0;
     thestatptr->itcount=0;
     thestatptr->totitcount=0;
     thestatptr->maxitcount=0;
     thestatptr->itzero=0;
     thestatptr->fcn=0;
      
     thestatptr->z = (val_type *)malloc(neq*(ns)*sizeof(val_type));
     thestatptr->li = (val_type *)malloc(neq*ns*sizeof(val_type));

     for (is=0; is<ns; is++)
          for (i=0; i<neq; i++)
               thestatptr->li[is*neq+i]=0.;
        
     return;

}


/************************************************************************************/
/* 					   					    */
/* NormalizedDistance: 	         	  					    */
/*                                        					    */
/*										    */
/************************************************************************************/

val_type NormalizedDistance (int neq,int ns,
                             toptions *options,val_type *z,val_type *zold)
{


/*---------------- Declarations -----------------------------------------*/

     int i,is,ix;
     val_type maxi,mayi,relerrors;
     val_type maxz,maxzold;

/* --------------- Implementation --------------------------------------*/
   	
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
               relerrors=FMAX(relerrors,FABS(z[ix]-zold[ix])/(mayi*options->rtol[i]+options->atol[i]));
          }

          maxi=FMAX(maxi,relerrors);

     }

     return maxi;

}

/************************************************************************************/
/* 					   					    */
/* RemoveDigitsFcn	         	  					    */
/*                                        					    */
/*										    */
/************************************************************************************/

void RemoveDigitsFcn (ode_sys *system,gauss_method *method,val_type *z, int m)
{

/* ---------- First initializations ------------------------------------*/

     int neq,ns;
    
     neq=system->neq;
     ns=method->ns;

/*---------------- Declarations -----------------------------------------*/

     int i,is,ix;
     val_type mY,mYY;

/* --------------- Implementation --------------------------------------*/

     for (i=0; i<neq; i++)
          for (is=0; is<ns; is++)
          { 
               ix=neq*is+i;
               mY=m*z[ix];
               mYY=mY+z[ix];
               z[ix]=mYY-mY;
           }

}

/************************************************************************************/
/* 					   					    */
/* Yi_init: 		         	  					    */
/*                                        					    */
/*										    */
/************************************************************************************/

int Yi_init(solution *u, val_type *z,ode_sys *system,
            gauss_method *method,solver_stat *thestatptr, toptions *options)
{


/* ---------- First initializations -------------------------------------*/

     int neq,ns;

     neq=system->neq;
     ns=method->ns;

/*------------- declarations -----------------------------------------*/ 

     int i,is,in,js,jsn;
     val_type *coef,*li;
     val_type sum;

/* ----------- implementation  --------------------------------------*/

 
     switch (options->approximation)
     {

     case 0:           

          for (is = 0; is<ns; is++)
               for (i = 0; i<neq; i++)
                    z[is*neq+i]=u->uu[i];       

     break;

     case 1:          
   
          li=thestatptr->li;
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
                }
          }  
     break;
          
     default:
          printf("approximation=%i error\n", options->approximation);
     break;

     }
 
     return(0);

}


/****************************************************************************************/
/* 											*/
/*   UpdateDmin  								        */
/* 											*/
/* 											*/
/****************************************************************************************/
void UpdateDMin (ode_sys *system,gauss_method *method,
                 bool *D0,bool *cont,val_type *DMin,val_type *Y,val_type *Yold)
{
/*----------------  First: initialization -------------------------------*/
     int neq,ns;

     neq=system->neq;
     ns=method->ns;

/*------ declarations --------------------------------------------------*/
     int i,is;
     val_type dY;

/* ----------- implementation  --------------------------------------*/ 

     *D0=true;
     *cont=false;

     for (is=0; is<ns; is++)
          for (i=0; i<neq; i++) 
          { 
               dY=FABS(Y[is*neq+i]-Yold[is*neq+i]);
               if (dY>0.)
               { 
                    *D0=false;
                    if (dY<DMin[is*neq+i])
                    {
                         DMin[is*neq+i]=dY;
                         *cont=true;
                    }
                   else
                    {
                         DMin[is*neq+i]=-1;
                    } 
               }
               else
               {
                    DMin[is*neq+i]=0;
               } 
          }
  
     return;

}


/****************************************************************************************/
/* 											*/
/*   Fixed_point_it  								        */
/* 											*/
/* 											*/
/****************************************************************************************/


void Fixed_point_it ( ode_sys *system, solution *u, val_type tn,val_type h, 
                      toptions *options,gauss_method *method,solver_stat *thestatptr)

{ 

/*----------------  First: initialization -------------------------------*/
     int neq,ns;

     neq=system->neq;
     ns=method->ns;

/*------ declarations --------------------------------------------------*/

     int i,is;
     bool D0,iter0;			
     val_type difftest,DMin[neq*ns];
     val_type *z;
     val_type zold[neq*ns];

/* ----------- implementation  --------------------------------------*/ 

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
           if (options->rdigits>0) RemoveDigitsFcn(system,method,z,options->mrdigits);
           UpdateDMin(system,method,&D0,&iter0,DMin,z,zold);
           thestatptr->itcount++;
     }

     if (thestatptr->itcount==MAXIT)
     { 
           printf("Break: step(MAXIT)=%i\n",thestatptr->itcount);
           thestatptr->convergence=FAIL; 
     }
 
     if (D0==false)
     {
           difftest=NormalizedDistance(neq,ns,options,z,zold);
           if (difftest>1.)
           {
                 printf("Lack of convegence of fixed-point iteration:\
                         step=%i,difftest=%lg\n",thestatptr->itcount,difftest);
                 thestatptr->convergence=FAIL;
           }
     }
     else
     {
           thestatptr->itzero++;
     }

     return;

}


/****************************************************************************************/
/*											*/
/*      It_Jacobi									*/
/*											*/
/****************************************************************************************/


int It_Jacobi   (ode_sys *system, solution *u, val_type tn,val_type h, 
                  gauss_method *method,solver_stat *thestatptr)
{

     int extern thread_count;

/* ---------- First initializations ------------------------------------*/

     int neq,ns;
     parameters *params;
     val_type *z,*li;
    
     neq=system->neq;
     ns=method->ns;
     params=&system->params;

     val_type fz[ns*neq];
     z=thestatptr->z;
     li=thestatptr->li;

/*------ declarations --------------------------------------------------*/

     int i,is,js,isn,in,jsn;
     val_type sum;
  
/* ----------- implementation  --------------------------------------*/ 


#ifdef PARALLEL
#      pragma omp parallel for num_threads(thread_count) private(isn)
#endif  
     for (is = 0; is<ns; is++)
     {
           isn=neq*is;
           params->ipar[0]=system->cod[0];
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


/****************************************************************************************/
/*											*/
/*      It_Jacobi_Classic								*/
/*											*/
/****************************************************************************************/


int It_Jacobi_Classic   (ode_sys *system, solution *u, val_type tn,val_type h, 
                         gauss_method *method,solver_stat *thestatptr)
{

     int extern thread_count;

/* ---------- First initializations ------------------------------------*/

     int neq,ns;
     parameters *params;
     val_type *z,*li;
    
     neq=system->neq;
     ns=method->ns;
     params=&system->params;

     val_type fz[ns*neq];
     z=thestatptr->z;
     li=thestatptr->li;

/*------ declarations --------------------------------------------------*/

     int i,is,js,isn,in,jsn;
     val_type sum;
  
/* ----------- implementation  --------------------------------------*/ 


#ifdef PARALLEL
#      pragma omp parallel for num_threads(thread_count) private(isn)
#endif  
     for (is = 0; is<ns; is++)
     {
           isn=neq*is;
           params->ipar[0]=system->cod[0];
           thestatptr->fcn++;
           system->f(neq,tn+method->hc[is],&z[isn],&fz[isn],params); 
           for (i=0; i<neq; i++) li[isn+i]=fz[isn+i]*method->hb[is];		
     }

     for (is = 0; is<ns; is++)
     { 
           for (i = 0; i<neq; i++)    
           {  
                 in=neq*is+i; 
//                 sum=method->m[ns*is]*li[i]+u->ee[i];    
                 sum=method->m[ns*is]*li[i];                      
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



/****************************************************************************************/
/*											*/
/*      It_Seidel									*/
/*											*/
/****************************************************************************************/


int It_Seidel   (ode_sys *system, solution *u, val_type tn,val_type h, 
                  gauss_method *method,solver_stat *thestatptr)
{

     int extern thread_count;

/* ---------- First initializations ------------------------------------*/

     int neq,ns;
     parameters *params;
     val_type *z,*li;
    
     neq=system->neq;
     ns=method->ns;
     params=&system->params;

     val_type fz[ns*neq];
     z=thestatptr->z;
     li=thestatptr->li;

/*------ declarations --------------------------------------------------*/

     int i,is,js,isn,in,jsn;
     val_type sum;
  
/* ----------- implementation  --------------------------------------*/ 


/* ----------- first part  -----------------------------------------*/ 
#ifdef PARALLEL
#      pragma omp parallel for num_threads(thread_count) private(isn)
#endif  
     for (is = 0; is<ns; is++)
     {
           isn=neq*is;
           params->ipar[0]=system->cod[0];
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


/* ----------- second part  --------------------------------------*/ 

#ifdef PARALLEL
#      pragma omp parallel for num_threads(thread_count) private(isn)
#endif  
     for (is = 0; is<ns; is++)
     {
           isn=neq*is;
           params->ipar[0]=system->cod[1];
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


 
void TheOutput(ode_sys *system,val_type t,solution *u,solver_stat *thestatptr,
               parameters *params,toptions *options,FILE *myfile)
{

/* ---------- First initializations ------------------------------------*/

     int neq;
     neq=system->neq;

#    if PREC ==2  //QUADRUPLEPRECISION
     int n;
     int width = 46;
     char buf[128];
#    endif

/*------ declarations --------------------------------------------------*/
     val_type DH;
     int i,node;
  
     struct rec
     {
     val_type t;
     val_type uu[neq];
     val_type ee[neq];
     };
  
     struct rec my_record;

/* ----------- implementation  -----------------------------------------*/ 

     node=neq;
     DH=(system->ham(node,u,params)-thestatptr->E0)/thestatptr->E0;
     if (fabs(DH)>thestatptr->MaxDE) thestatptr->MaxDE=FABS(DH);

     if (((thestatptr->stepcount % options->sampling) == 0) || (thestatptr->laststep))
     {
           thestatptr->nout++;
#    if PREC ==2  //QUADRUPLEPRECISION
           printf("Energy Quadruple, %i,%i,",thestatptr->stepcount,thestatptr->itcount);
           n = quadmath_snprintf(buf, sizeof buf, "%+-#*.30Qe", width, DH);
           if ((size_t) n < sizeof buf) printf("%s\n",buf);
#    else  //DOUBLEPRECISION
           printf("%i, %lg\n",thestatptr->stepcount,DH);
#    endif 


           my_record.t=t;
           for (i=0; i<neq;i++) 
           {
                 my_record.uu[i]=u->uu[i];
                 my_record.ee[i]=u->ee[i];
           }
           fwrite(&my_record, sizeof(struct rec), 1, myfile);
  
     } 

  return;

}

void TheOutput2(ode_sys *system,val_type t,solution *u,solution *u2,solver_stat *thestatptr,
                    parameters *params,toptions *options,FILE *myfile)
{

/* ---------- First initializations ------------------------------------*/

     int neq;
     neq=system->neq;

/*------ declarations --------------------------------------------------*/

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

#    if PREC ==2  //QUADRUPLEPRECISION
     int n;
     int width = 46;
     char buf[128];
#    endif

/* ----------- implementation  -----------------------------------------*/ 

     node=neq;
     DH=(system->ham(node,u,params)-thestatptr->E0)/thestatptr->E0;
     if (fabs(DH)>thestatptr->MaxDE) thestatptr->MaxDE=fabs(DH);

     if (((thestatptr->stepcount % options->sampling) == 0) || (thestatptr->laststep))
     {
           thestatptr->nout++;
#     if PREC ==2  //QUADRUPLEPRECISION
           printf("Energy Quadruple, %i,",thestatptr->stepcount);
           n = quadmath_snprintf(buf, sizeof buf, "%+-#*.30Qe", width, DH);
           if ((size_t) n < sizeof buf) printf("%s\n",buf);
#     else     //DOUBLEPRECISION
           printf("%i, %lg\n",thestatptr->stepcount,DH);
#     endif 

     my_record.t=t;
     for (i=0; i<neq;i++) 
     {
           my_record.uu[i]=u->uu[i];
           my_record.ee[i]=u->ee[i];
           my_record.est[i]=(u->uu[i]-u2->uu[i])+(u->ee[i]-u2->ee[i]);}
           fwrite(&my_record, sizeof(struct rec), 1, myfile);  
     } 

     return;

}



/****************************************************************************************/
/* 											*/
/*   RGK:  									        */
/* 											*/
/* 											*/
/****************************************************************************************/

void RKG (gauss_method *gsmethod, solution *u,ode_sys *system,toptions *options,
          void RKG_Step (), solver_stat *thestatptr)

{

/* ---------- First initializations ------------------------------------*/

     int neq,ns;
     parameters *params;
     
     neq=system->neq;
     ns=gsmethod->ns;
     params=&system->params;

/*------ declarations --------------------------------------------------*/

     FILE *myfile;
  
     int i,is,isn,ix;
     int istep,nstep;
     val_type tn;             
     val_type u0,sum,aux;  
     val_type *z,*li;

     z=thestatptr->z;
     li=thestatptr->li;


#    if PREC ==2  //QUADRUPLEPRECISION
     int n;
     int width = 46;
     char buf[128];
#    endif

/* ----------- implementation  ---------------------------------------*/ 

     myfile = fopen(thestatptr->filename,"wb"); 

/* ----------- initial energi (E0)   ---------------------------------*/

     thestatptr->E0=system->ham(neq,u,params);

#if PREC ==2  //QUADRUPLEPRECISION
     printf("Initial energy:");
     n = quadmath_snprintf(buf, sizeof buf, "%+-#*.30Qe", width, thestatptr->E0);
     if ((size_t) n < sizeof buf) printf("%s\n",buf);

#else     //DOUBLEPRECISION
     printf("Initial energy=%lg\n", thestatptr->E0);
#endif 

     thestatptr->MaxDE=0.; 
     tn=options->t0;
     nstep=((options->t1)-(options->t0))/(options->h);
     if (nstep*(options->h)+(options->t0)<options->t1)  {nstep++;}    
//     nstep=1;options->t1=nstep*(options->h);
  
     thestatptr->convergence=SUCCESS;

     thestatptr->nout=0;
 
#ifdef IOUT
     TheOutput(system,tn,u,thestatptr,params,options,myfile);
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

     /* ----------------yn+1=yn+(Sum Li+e_n) --------------------------------------*/

          for (i = 0; i<neq; i++)
          { 
               ix=0;
               is=gsmethod->orderedindices[ix];
               sum=li[neq*is+i]+u->ee[i];			

               for (ix =1 ; ix<ns; ix++)
               {
                    is=gsmethod->orderedindices[ix];
                    isn=neq*is+i;
                    sum+=li[isn];
               } 

     /* ----------------compensated sumation ------------------------------------*/  

               u0=u->uu[i];
               u->uu[i]=u0+sum;
               aux=u->uu[i]-u0;
               u->ee[i]=sum-aux;
      	       	       
          }

          tn=(istep+1)*(options->h);

          thestatptr->stepcount++;		
          (thestatptr->totitcount)+=(thestatptr->itcount); 
          if ((thestatptr->itcount)>(thestatptr->maxitcount))
               (thestatptr->maxitcount)=(thestatptr->itcount);

#ifdef IOUT
          TheOutput(system,tn,u,thestatptr,params,options,myfile);
#endif 

          if ((tn+options->h)>=options->t1)
          { 
               options->h=options->t1-tn;
               thestatptr->laststep=true;
          } 


     }

     fclose(myfile);

     return;

}

/****************************************************************************************/
/* 											*/
/*   RKG2:  double integration (sequentially).  				        */
/*        									        */
/* 											*/
/****************************************************************************************/

void RKG2 (gauss_method *gsmethod, gauss_method *gsmethod2,solution *u,solution *u2,
           ode_sys *system,toptions *options,toptions *options2,
           void RKG_Step (), solver_stat *thestatptr,solver_stat *thestatptr2)
 
{

/* ---------- First initializations ------------------------------------*/

     int neq,ns;
     parameters *params;

     neq=system->neq;
     ns=gsmethod->ns;
     params=&system->params;

/*------ declarations --------------------------------------------------*/

     FILE *myfile;
  
     int i,is,isn,ix,ksw;
     int istep,nstep;
     val_type tn;             
     val_type sum,aux;
     val_type u0[neq];     
     val_type *z,*li,*z2,*li2;
     bool initwithfirstintegration;         

     z=thestatptr->z;
     li=thestatptr->li;
     z2=thestatptr2->z;
     li2=thestatptr2->li;

#    if PREC ==2  //QUADRUPLEPRECISION
     int n;
     int width = 46;
     char buf[128];
#    endif

/* ----------- implementation  ---------------------------------------*/ 

     myfile = fopen(thestatptr->filename,"wb"); 

/* ----------- initial energi (E0)   ---------------------------------*/

     thestatptr->E0=system->ham(neq,u,params);

#if PREC ==2  //QUADRUPLEPRECISION
     printf("Initial energy:");
     n = quadmath_snprintf(buf, sizeof buf, "%+-#*.30Qe", width, thestatptr->E0);
     if ((size_t) n < sizeof buf) printf("%s\n",buf);

#else     //DOUBLEPRECISION
     printf("Initial energy=%lg\n", thestatptr->E0);
#endif 

     thestatptr->MaxDE=0.; 
  
     tn=options->t0;
     nstep=((options->t1)-(options->t0))/(options->h);
     if (nstep*(options->h)+(options->t0)<options->t1) 
           nstep++;    
//  nstep=10;options->t1=nstep*(options->h);
  
     thestatptr->convergence=SUCCESS;
     initwithfirstintegration=true;
     ksw=0;

     thestatptr->nout=0;

#ifdef IOUT
     TheOutput2(system,tn,u,u2,thestatptr,params,options,myfile);
#endif 

     for(istep=0; istep<nstep; istep++) 
     {     

          /* ---------------- First integration -----------------------------------------*/

          Yi_init (u,z,system,gsmethod,thestatptr,options);   
          RKG_Step (system,u,tn,options->h,options,gsmethod,thestatptr);
          if (thestatptr->convergence==FAIL)
          {
               printf("Stop Fail. step=%i\n", istep);
               nstep=istep;
          } 	      

          for (i = 0; i<neq; i++)
          { 
               ix=0;
               is=gsmethod->orderedindices[ix];
               sum=li[neq*(is)+i]+u->ee[i];			

               for (ix = 1 ; ix<ns; ix++)
               {
                    is=gsmethod->orderedindices[ix];
                    isn=neq*is+i;
                    sum+=li[isn];
               } 

               u0[i]=u->uu[i];
               u->uu[i]=u0[i]+sum;
               aux=u->uu[i]-u0[i];
               u->ee[i]=sum-aux;
      	       	       
          }

          /* ---------------- Second integration -----------------------------------------*/

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
          if (thestatptr2->convergence==FAIL)
          {
               printf("Stop Fail. step=%i\n", istep);
               nstep=istep;
          } 	
  
          for (i = 0; i<neq; i++)
          { 
               ix=0;
               is=gsmethod2->orderedindices[ix];
               sum=li2[neq*(is)+i]+u2->ee[i];			

               for (ix = 1 ; ix<ns; ix++)
               {
                    is=gsmethod2->orderedindices[ix];
                    isn=neq*is+i;
                    sum+=li2[isn];
               } 

               u0[i]=u2->uu[i];
               u2->uu[i]=u0[i]+sum;
               aux=u2->uu[i]-u0[i];
               u2->ee[i]=sum-aux;
      	       	       
          }
   
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
          TheOutput2(system,tn,u,u2,thestatptr,params,options,myfile);
#endif 

          if ((tn+options->h)>=options->t1)
          { 
               options->h=options->t1-tn;
               thestatptr->laststep=true;
          } 

      }


      fclose(myfile);

      return;

}


/****************************************************************************************/
/* 											*/
/*   select_gauss				 				        */
/*        									        */
/* 											*/
/****************************************************************************************/

void select_gauss(gauss_method *gsmethodptr, gauss_method *gsmethod2ptr,solution *uptr,solution *u2ptr,
           ode_sys *systemptr,toptions *optionsptr,toptions *options2ptr,
           solver_stat *thestatptr,solver_stat *thestat2ptr)
{
    int i;

    switch (optionsptr->algorithm)
    {  
    case  1:  /* Jacobi */
          optionsptr->iteration[0]=It_Jacobi;
          systemptr->cod[0]=0;
          RKG (gsmethodptr,uptr,systemptr,optionsptr,Fixed_point_it,thestatptr);                             
    break;  

    case  9:  /* Jacobi_Classic */
          optionsptr->iteration[0]=It_Jacobi_Classic;
          systemptr->cod[0]=0;
          RKG (gsmethodptr,uptr,systemptr,optionsptr,Fixed_point_it,thestatptr);                             
    break; 

    case  11:  /* Seidel */
          optionsptr->iteration[0]=It_Seidel;
          systemptr->cod[0]=1;
          systemptr->cod[1]=2;
	  RKG (gsmethodptr,uptr,systemptr,optionsptr,Fixed_point_it,thestatptr);
    break;  

    case 21: /* RKG2: sequential execution (jacobi) */
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

          GaussCoefficients(gsmethod2ptr,options2ptr);
 	  InitStat(systemptr,gsmethod2ptr,thestat2ptr);

          RKG2 (gsmethodptr,gsmethod2ptr,uptr,u2ptr,systemptr,optionsptr,options2ptr,
                Fixed_point_it,thestatptr,thestat2ptr);  

    break;

    case 22: /* RKG2: sequential execution (seidel) */
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

          optionsptr->iteration[0]=It_Seidel;
          options2ptr->iteration[0]=It_Seidel;
          systemptr->cod[0]=1;
          systemptr->cod[1]=2;   
          gsmethod2ptr->ns=gsmethodptr->ns;

          GaussCoefficients(gsmethod2ptr,options2ptr);
 	  InitStat(systemptr,gsmethod2ptr,thestat2ptr);

          RKG2 (gsmethodptr,gsmethod2ptr,uptr,u2ptr,systemptr,optionsptr,options2ptr,
                Fixed_point_it,thestatptr,thestat2ptr);  

    break;
                                    
    default:
          printf("Error: incorrect algorithm\n");
    break;

         
    } 
}


