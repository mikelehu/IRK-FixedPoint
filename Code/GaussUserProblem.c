/*------------------------------------------------------------------------------*/
/*										*/
/*    GaussUserProblem.c							*/
/*	Functions:								*/
/*	   Ode1()= OdePendulum():						*/
/*	   Ham1()= HamPendulum():						*/
/*         Ode2()= OdeNbody():							*/
/*	   Ham2()= HamNBody():							*/
/*										*/
/*	   Ode3,...,Ode10 : empty functions					*/
/*	   Ham3,...,Ham10 : empty functions					*/
/*										*/
/*         Ode11= OdePendulumX(): ideal integrator				*/
/*         Ode12= OdeNBodyX():    ideal integrator				*/
/*										*/
/*------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <def.h>
#include <quadmath.h>


/*------------------------------------------------------------------------------*/
/*										*/
/*       OdePendulum Problem: 							*/
/*	    Ode1()=OdePendulum							*/
/*	    Ham1()=HamPendulum():						*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode1 (int neq, val_type t,val_type *u,val_type *f,parameters *params)

{
    
 /* ---------- First initializations -----------------------------------------*/

 /*------ declarations -------------------------------------------------*/ 

     val_type C1,C2,C3,C4,C5,C6,C7;
     val_type Q1,Q2,P1,P2,Q12;
     val_type cosQ1Q2,sinQ1Q2,sinQ1,sinQ2;
     val_type aux0,aux1,dQ1aux1,aux2,aux3,aux4;

 /* ----------- implementation  ---------------------------------------*/

     C1=params->rpar[0];
     C2=params->rpar[1];
     C3=params->rpar[2];
     C4=params->rpar[3];
     C5=params->rpar[4];
     C6=params->rpar[5];
     C7=params->rpar[6];

     Q1=u[0];
     Q2=u[1];
     P1=u[2];
     P2=u[3];
     Q12=(u[0]-u[1]);
           
     cosQ1Q2= COS(Q12);
     sinQ1Q2= SIN(Q12);
     sinQ1= SIN(Q1);
     sinQ2= SIN(Q2);

     aux0= C5*sinQ1Q2*sinQ1Q2;
     aux1= C4+aux0;
     dQ1aux1= 2*C5*sinQ1Q2*cosQ1Q2;
     aux2= aux1*aux1;
     aux3= C3*cosQ1Q2;
     aux4= (-1/aux2)*(C1*P1*P1+C2*P2*P2+P1*P2*aux3)*dQ1aux1-(C3*P1*P2*sinQ1Q2)/aux1;

     f[0]= (2*C1*P1+aux3*P2)/aux1;
     f[1]= (2*C2*P2+aux3*P1)/aux1;
     f[2]=-(aux4+C6*sinQ1);
     f[3]=-(-aux4+C7*sinQ2);               
                        
     return ;

}



val_type Ham1 (int neq,solution *u,parameters *params)
{

/* ---------- First initializations -----------------------------------------*/

     val_type *uu;
     uu=u->uu;

/*------ declarations -------------------------------------------------*/ 

     val_type C1,C2,C3,C4,C5,C6,C7;
     val_type Q1,Q2,P1,P2,Q12;
     val_type cosQ1Q2,sinQ1Q2,cosQ1,cosQ2;
     val_type H;

/* ----------- implementation  ---------------------------------------*/

     C1=params->rpar[0];
     C2=params->rpar[1];
     C3=params->rpar[2];
     C4=params->rpar[3];
     C5=params->rpar[4];
     C6=params->rpar[5];
     C7=params->rpar[6];

     Q1=uu[0];
     Q2=uu[1];
     P1=uu[2];
     P2=uu[3];
     Q12=Q1-Q2;
           
     cosQ1Q2= COS(Q12);
     sinQ1Q2= SIN(Q12);
     cosQ1= COS(Q1);
     cosQ2= COS(Q2);

     H= ((C1*P1*P1+C2*P2*P2+C3*P1*P2*cosQ1Q2) / (C4+C5*sinQ1Q2*sinQ1Q2))
        - C6*cosQ1-C7*cosQ2;

     return(H);

}



/*------------------------------------------------------------------------------*/
/*										*/
/*       NBody Problem: 							*/
/*	    Ode2()=OdeNBody							*/
/*	    Ham2()=HamNBody():							*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode2 (int neq, val_type t,val_type *u,val_type *f,parameters *params)

{
/* ---------- First initializations ------------------------------------*/
  
     int dim;
     dim=3;
    
/*------ declarations -------------------------------------------------*/ 

     int i,id,i1,i2,j,j1,j2;
     int nd,nbody;
     val_type d3,qij;
     val_type *Gm;

/* ----------- implementation  ---------------------------------------*/

     nbody=neq/(2*dim);
     nd=neq/2;
     Gm=params->rpar;
    
     switch (params->ipar[0])
     {
     case 1: /*OdePlanetsq*/

           for (i=0; i<nbody; i++)
           {
                 i1=i*dim;
                 i2=nd+i1;
                 for (id=0; id<dim; id++) f[i1+id]=u[i2+id];
           }
              
     break;

     case 2: /*OdePlanetsv*/      

           for (i=0; i<nbody; i++)
           {
                 i1=i*dim;
                 i2=nd+i1;
                 for (id=0; id<dim; id++) f[i2+id]=0.;
           }
        
           for (i=0; i<nbody; i++)
           {
                 i1=i*dim;
                 i2=nd+i1;

                 for (j=i+1; j<nbody; j++)
                 {
                       j1=j*dim;
                       j2=nd+j1;
                       d3=0.;
                       for (id=0; id<dim; id++) d3+=POW(u[i1+id]-u[j1+id],2);
                       d3=POW(d3,3./2);

                       for (id=0; id<dim; id++) 
                       {
                             qij=(u[i1+id]-u[j1+id])/d3;
                             f[i2+id]-=Gm[j]*qij;
                             f[j2+id]+=Gm[i]*qij;
                       }   
                 }
           }
         
     break;
                                     
     default:     
        
           for (i=0; i<nbody; i++)
           {
                 i1=i*dim;
                 i2=nd+i1;
                 for (id=0; id<dim; id++) 
                 {
                       f[i1+id]=u[i2+id];
                       f[i2+id]=0.;
                 }
           }

           for (i=0; i<nbody; i++)
           {
                 i1=i*dim;
                 i2=nd+i1;
                 for (j=i+1; j<nbody; j++)
                 {
                       j1=j*dim;
                       j2=nd+j1;
                       d3=0.;
                       for (id=0; id<dim; id++) d3+=POW(u[i1+id]-u[j1+id],2);
                       d3=POW(d3,3./2);

                       for (id=0; id<dim; id++) 
                       {
                             qij=(u[i1+id]-u[j1+id])/d3;
                             f[i2+id]-=Gm[j]*qij;
                             f[j2+id]+=Gm[i]*qij;
                       }   
                 }

           }
 
     break;

}


    return ;

}


val_type Ham2 (int neq,solution *u,parameters *params)
{

/* ---------- First initializations ------------------------------------*/
  
     int dim;
     dim=3;

/*------ declarations -------------------------------------------------*/ 

     val_type *uu;  
     int i,j,i1,j1;
     int nbody,node,nd;
     val_type *d,*Gm;
     val_type H,Pot;
 
/* ----------- implementation  ----------------------------------------*/
     uu=u->uu;
     Gm=params->rpar;
     nbody=neq/(2*dim);
     d = malloc(nbody*nbody*sizeof(val_type));
     node=neq;
     nd=node/2;

     H=0.;
 
     for (i=0; i<nbody; i++)
     {
           i1=dim*i;
           for (j=i+1; j<nbody; j++)
           {
                 j1=dim*j;
                 d[i*nbody+j]=SQRT(POW(uu[i1]-uu[j1],2)+
                                   POW(uu[i1+1]-uu[j1+1],2)+POW(uu[i1+2]-uu[j1+2],2));
                 d[j*nbody+i]=d[i*nbody+j];
           }
     }

     for (i=0; i<nbody; i++)
     {
           i1=nd+dim*i;
           H+=Gm[i]*(POW(uu[i1],2)+POW(uu[i1+1],2)+POW(uu[i1+2],2));
     }

     H=H/2.;
     Pot=0;

     for (i=0; i<nbody-1; i++)
          for (j=i+1; j<nbody; j++)
               Pot+=Gm[i]*Gm[j]/d[i*nbody+j];


     H=H-Pot;
 
     free(d);
     return(H);

}





/*------------------------------------------------------------------------------*/
/*										*/
/*        Problem-3:	 							*/
/*	    Ode3()=								*/
/*	    Ham3()=								*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode3 (int neq, val_type t,val_type *u,val_type *f,parameters *params)
{
     return;
}

val_type Ham3 (int neq,solution *u,parameters *params)
{
     return(0.);

}


/*------------------------------------------------------------------------------*/
/*										*/
/*        Problem-4:	 							*/
/*	    Ode4()=								*/
/*	    Ham4()=								*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode4 (int neq, val_type t,val_type *u,val_type *f,parameters *params)
{
     return;
}

val_type Ham4 (int neq,solution *u,parameters *params)
{
     return(0.);

}

/*------------------------------------------------------------------------------*/
/*										*/
/*        Problem-5:	 							*/
/*	    Ode5()=								*/
/*	    Ham5()=								*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode5 (int neq, val_type t,val_type *u,val_type *f,parameters *params)
{
     return;
}

val_type Ham5 (int neq,solution *u,parameters *params)
{
     return(0.);

}


/*------------------------------------------------------------------------------*/
/*										*/
/*        Problem-6:	 							*/
/*	    Ode6()=								*/
/*	    Ham6()=								*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode6 (int neq, val_type t,val_type *u,val_type *f,parameters *params)
{
     return;
}

val_type Ham6 (int neq,solution *u,parameters *params)
{
     return(0.);

}


/*------------------------------------------------------------------------------*/
/*										*/
/*        Problem-7:	 							*/
/*	    Ode7()=								*/
/*	    Ham7()=								*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode7 (int neq, val_type t,val_type *u,val_type *f,parameters *params)
{
     return;
}

val_type Ham7 (int neq,solution *u,parameters *params)
{
     return(0.);

}


/*------------------------------------------------------------------------------*/
/*										*/
/*        Problem-8:	 							*/
/*	    Ode8()=								*/
/*	    Ham8()=								*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode8 (int neq, val_type t,val_type *u,val_type *f,parameters *params)
{
     return;
}

val_type Ham8 (int neq,solution *u,parameters *params)
{
     return(0.);

}



/*------------------------------------------------------------------------------*/
/*										*/
/*        Problem-9:	 							*/
/*	    Ode9()=								*/
/*	    Ham9()=								*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode9 (int neq, val_type t,val_type *u,val_type *f,parameters *params)
{
     return;
}

val_type Ham9 (int neq,solution *u,parameters *params)
{
     return(0.);

}



/*------------------------------------------------------------------------------*/
/*										*/
/*        Problem-10:	 							*/
/*	    Ode10()=								*/
/*	    Ham10()=								*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode10 (int neq, val_type t,val_type *u,val_type *f,parameters *params)
{
     return;
}

val_type Ham10 (int neq,solution *u,parameters *params)
{
     return(0.);

}


/*------------------------------------------------------------------------------*/
/*										*/
/*       OdePendulum Problem: 							*/
/*		  Ode11= OdePendulumX (Ideal integrator)			*/
/*										*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode11 (int neq, val_type t,val_type *u,val_type *f,parameters *params)

{
 /* --------------------------------------------------------------------------*/
 /* specific Odependulum() for "ideal integrator" esperiment                  */
 /*     - input:      quadruple---> double                                    */
 /*     - operations: double                                                  */
 /*     - output:     double ---> quadruple.                                  */
 /* --------------------------------------------------------------------------*/

 /* ---------- First initializations -----------------------------------------*/

 /*------ declarations -------------------------------------------------*/ 

     int i;

     double fdouble[neq];

     double C1,C2,C3,C4,C5,C6,C7;
     double Q1,Q2,P1,P2,Q12;
     double cosQ1Q2,sinQ1Q2,sinQ1,sinQ2;
     double aux0,aux1,dQ1aux1,aux2,aux3,aux4;

 /* ----------- implementation  ---------------------------------------*/

     C1=params->rpar[0];
     C2=params->rpar[1];
     C3=params->rpar[2];
     C4=params->rpar[3];
     C5=params->rpar[4];
     C6=params->rpar[5];
     C7=params->rpar[6];

     Q1=u[0];
     Q2=u[1];
     P1=u[2];
     P2=u[3];
     Q12=(u[0]-u[1]);
            
     cosQ1Q2= cos(Q12);
     sinQ1Q2= sin(Q12);
     sinQ1= sin(Q1);
     sinQ2= sin(Q2);

     aux0= C5*sinQ1Q2*sinQ1Q2;
     aux1= C4+aux0;
     dQ1aux1= 2*C5*sinQ1Q2*cosQ1Q2;
     aux2= aux1*aux1;
     aux3= C3*cosQ1Q2;
     aux4= (-1/aux2)*(C1*P1*P1+C2*P2*P2+P1*P2*aux3)*dQ1aux1-(C3*P1*P2*sinQ1Q2)/aux1;

     fdouble[0]= (2*C1*P1+aux3*P2)/aux1;
     fdouble[1]= (2*C2*P2+aux3*P1)/aux1;
     fdouble[2]=-(aux4+C6*sinQ1);
     fdouble[3]=-(-aux4+C7*sinQ2);               
              
     for (i=0; i<neq; i++) f[i]=fdouble[i];
          
     return ;

}


/*------------------------------------------------------------------------------*/
/*										*/
/*       NBody Problem:								*/
/*	       Ode12=OdeNBodyX (Ideal integrator)				*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode12 (int neq, val_type t,val_type *u,val_type *f,parameters *params)

{

 /* --------------------------------------------------------------------------*/
 /* specific Odependulum() for "ideal integrator" esperiment                  */
 /*     - input:      quadruple---> double                                    */
 /*     - operations: double                                                  */
 /*     - output:     double ---> quadruple.                                  */
 /* --------------------------------------------------------------------------*/


/* ---------- First initializations ------------------------------------*/
  
     int dim,nd,nbody;
     dim=3;
     nbody=neq/(2*dim);
     nd=neq/2;
    
/*------ declarations -------------------------------------------------*/ 

     int i,id,i1,i2,j,j1,j2;
     double d3,qij;
     val_type *Gm;

     double fdouble[neq],Udouble[neq],Gmdouble[nbody];

/* ----------- implementation  ---------------------------------------*/

     Gm=params->rpar;

     for (i=0;i<nbody;i++) Gmdouble[i]=Gm[i];
     for (i=0;i<neq;i++) Udouble[i]=u[i];
    
     switch (params->ipar[0])
     {
     case 1: /*OdePlanetsq*/

           for (i=0; i<nbody; i++)
           {
                 i1=i*dim;
                 i2=nd+i1;
                 for (id=0; id<dim; id++) fdouble[i1+id]=Udouble[i2+id];
           }
              
     break;

     case 2: /*OdePlanetsv*/      

           for (i=0; i<nbody; i++)
           {
                 i1=i*dim;
                 i2=nd+i1;
                 for (id=0; id<dim; id++) fdouble[i2+id]=0.;
           }
        
           for (i=0; i<nbody; i++)
           {
                 i1=i*dim;
                 i2=nd+i1;

                 for (j=i+1; j<nbody; j++)
                 {
                       j1=j*dim;
                       j2=nd+j1;
                       d3=0.;
                       for (id=0; id<dim; id++) d3+=POW(Udouble[i1+id]-Udouble[j1+id],2);
                       d3=POW(d3,3./2);

                       for (id=0; id<dim; id++) 
                       {
                             qij=(Udouble[i1+id]-Udouble[j1+id])/d3;
                             fdouble[i2+id]-=Gmdouble[j]*qij;
                             fdouble[j2+id]+=Gmdouble[i]*qij;
                       }   
                 }
           }
         
     break;
                                     
     default:     
        
           for (i=0; i<nbody; i++)
           {
                 i1=i*dim;
                 i2=nd+i1;
                 for (id=0; id<dim; id++) 
                 {
                       fdouble[i1+id]=Udouble[i2+id];
                       fdouble[i2+id]=0.;
                 }
           }

           for (i=0; i<nbody; i++)
           {
                 i1=i*dim;
                 i2=nd+i1;
                 for (j=i+1; j<nbody; j++)
                 {
                       j1=j*dim;
                       j2=nd+j1;
                       d3=0.;
                       for (id=0; id<dim; id++) d3+=POW(Udouble[i1+id]-Udouble[j1+id],2);
                       d3=POW(d3,3./2);

                       for (id=0; id<dim; id++) 
                       {
                             qij=(Udouble[i1+id]-Udouble[j1+id])/d3;
                             fdouble[i2+id]-=Gmdouble[j]*qij;
                             fdouble[j2+id]+=Gmdouble[i]*qij;
                       }   
                 }

           }
 
     break;

}

    for (i=0; i<neq; i++) f[i]=fdouble[i];

    return ;

}




