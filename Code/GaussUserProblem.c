/*------------------------------------------------------------------------------*/
/*										*/
/*    GaussUserProblem.c							*/
/*	Functions:								*/
/*	   Ode1()= OdePendulum():						*/
/*	   Ham1()= HamPendulum():						*/
/*         Ode2()= OdeNbody():							*/
/*	   Ham2()= HamNBody():							*/
/*	   Ode3()= empty:				         		*/
/*	   Ham3()= empty:							*/
/*	   Ode4()= OdePendulumStiff():	Gaizki !!!				*/
/*	   Ham4()= HamPendulumStiff():	Gaizki !!!				*/
/*	   Ode5()= OdePendulumV2():  NEW 13.12.2016				*/
/*	   Ham5()= HamPendulumV2():  NEW 13.12.2016				*/
/*										*/
/*	   Ode6,...,Ode10 : empty functions					*/
/*	   Ham6,...,Ham10 : empty functions					*/
/*										*/
/*         Ode11= OdePendulumX(): ideal integrator				*/
/*         Ode12= OdeNBodyX():    ideal integrator                              */
/*         Ode14= OdePendulumStiffX()    ideal integrator  			*/
/*         Ode15= OdePendulumX(): ideal integrator	NEW 13.12.2016		*/
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

void Ode1 (const int neq,const val_type t,const val_type *u,
           val_type *f, const parameters *params)

{
    
/* ---------- First initializations -----------------------------------------*/

/*------ declarations -------------------------------------------------------*/ 

     val_type C1,C2,C3,C4,C5,C6,C7;
     val_type Q1,Q2,P1,P2,Q12;
     val_type cosQ1Q2,sinQ1Q2,sinQ1,sinQ2;
     val_type aux0,aux1,dQ1aux1,aux2,aux3,aux4;

/* ----------- implementation  ----------------------------------------------*/

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

     Q12=(Q1-Q2);           
     cosQ1Q2= COS(Q12);
     sinQ1Q2= SIN(Q12);
     sinQ1= SIN(Q1);
     sinQ2= SIN(Q2);

     aux0= C5*sinQ1Q2*sinQ1Q2;
     aux1= (C4+aux0);
     dQ1aux1= 2*C5*sinQ1Q2*cosQ1Q2;
     aux2= aux1*aux1;
     aux3= C3*cosQ1Q2;

     aux4= (-1/aux2)*(C1*P1*P1+C2*P2*P2+P1*P2*aux3)*dQ1aux1-(C3*P1*P2*sinQ1Q2)/aux1;

     f[0]= (2*C1*P1+aux3*P2)/aux1;
     f[1]= (2*C2*P2+aux3*P1)/aux1;
     f[2]=-(aux4+C6*sinQ1); 
     f[3]=(aux4-C7*sinQ2); 

/*     printf("\nOdefun1\n");
     printf("%.20lg,%.20lg,%.20lg,%.20lg\n", f[0],f[1],f[2],f[3]);    
*/
                      
     return ;

}



val_type Ham1 (const int neq, const solution *u, const parameters *params)
{

/* ---------- First initializations -----------------------------------------*/

     val_type *uu;
     uu=u->uu;

/*------ declarations -------------------------------------------------------*/ 

     val_type C1,C2,C3,C4,C5,C6,C7;
     val_type Q1,Q2,P1,P2,Q12;
     val_type cosQ1Q2,sinQ1Q2,cosQ1,cosQ2;
     val_type H;

/* ----------- implementation  ----------------------------------------------*/

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

void Ode2 (const int neq,const val_type t,const val_type *u,
           val_type *f, const parameters *params)

{
/* ---------- First initializations -----------------------------------------*/
  
     int dim;
     dim=3;
    
/*------ declarations -------------------------------------------------------*/ 

     int i,id,i1,i2,j,j1,j2;
     int ix,jx;
     int nd,nbody;
     val_type d3,qij;
     val_type *Gm;

/* ----------- implementation  ----------------------------------------------*/

/*   params->ipar=orderplanet (ascending order)*/

     nbody=neq/(2*dim);
     nd=neq/2;
     Gm=params->rpar;
    
     switch (params->eval)
     {
     case 1: /*OdePlanetsq*/

           for (ix=0; ix<nbody; ix++)
           {
                 i=params->ipar[ix];		 
                 i1=i*dim;
                 i2=nd+i1;
                 for (id=0; id<dim; id++) f[i1+id]=u[i2+id];
           }
              
     break;

     case 2: /*OdePlanetsv*/      

           for (ix=0; ix<nbody; ix++)
           {
                 i=params->ipar[ix];
                 i1=i*dim;
                 i2=nd+i1;
                 for (id=0; id<dim; id++) f[i2+id]=0.;
           }
        
           for (ix=0; ix<nbody; ix++)
           {
                 i=params->ipar[ix];
                 i1=i*dim;
                 i2=nd+i1;

                 for (jx=ix+1; jx<nbody; jx++)
                 {
                       j=params->ipar[jx];
                       j1=j*dim;
                       j2=nd+j1;
                       d3=0.;
                       for (id=0; id<dim; id++) d3+=(u[i1+id]-u[j1+id])*(u[i1+id]-u[j1+id]);
                       d3=SQRT(d3)*d3;
//                       for (id=0; id<dim; id++) d3+=POW(u[i1+id]-u[j1+id],2);
//                       d3=POW(d3,3./2);

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
        
           for (ix=0; ix<nbody; ix++)
           {
                 i=params->ipar[ix];
//                 i=ix;
                 i1=i*dim;
                 i2=nd+i1;
                 for (id=0; id<dim; id++) 
                 {
                       f[i1+id]=u[i2+id];
                       f[i2+id]=0.;
                 }
           }

           for (ix=0; ix<nbody; ix++)
           {
                 i=params->ipar[ix];
                 i1=i*dim;
                 i2=nd+i1;
                 for (jx=ix+1; jx<nbody; jx++)
                 {
                       j=params->ipar[jx];
                       j1=j*dim;
                       j2=nd+j1;
                       d3=0.;
                       for (id=0; id<dim; id++) d3+=(u[i1+id]-u[j1+id])*(u[i1+id]-u[j1+id]);
                       d3=SQRT(d3)*d3;
//                       for (id=0; id<dim; id++) d3+=POW(u[i1+id]-u[j1+id],2);
//                       d3=POW(d3,3./2);

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

/*
    printf("\nOdefun2\n");
    for (i=0; i<neq;i++) printf("%lg,", f[i]);
    printf("\n");
*/

    return ;

}



val_type Ham2 (const int neq, const solution *u, const parameters *params)
{

/* ---------- First initializations -----------------------------------------*/
  
     int dim;
     dim=3;

/*------ declarations -------------------------------------------------------*/ 

     val_type *uu;  
     int i,j,i1,j1;
     int nbody,node,nd;
     val_type *d,*Gm;
     val_type H,Pot;
 
/* ----------- implementation  ----------------------------------------------*/
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
/*	    Ode3()=:								*/
/*	    Ham3()=:								*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode3 (const int neq,const val_type t,const val_type *u,
           val_type *f, const parameters *params)
{
 
                        
     return ;

}

val_type Ham3 (const int neq, const solution *u, const parameters *params)
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

void Ode4 (const int neq,const val_type t,const val_type *u,
           val_type *f, const parameters *params)

{

 /* ---------- First initializations -----------------------------------------*/

 /*------ declarations -------------------------------------------------*/ 

     val_type C1,C2,C3,C4,C5,C6,C7,k;
     val_type Q1,Q2,P1,P2;
     val_type cosQ2,cos2Q2;
     val_type aux1,aux2,aux3,aux4,aux5,aux6,aux7;

 /* ----------- implementation  ---------------------------------------*/

     C1=params->rpar[0];
     C2=params->rpar[1];
     C3=params->rpar[2];
     C4=params->rpar[3];
     C5=params->rpar[4];
     C6=params->rpar[5];
     C7=params->rpar[6];
     k=params->rpar[7];

     Q1=u[0];
     Q2=u[1];
     P1=u[2];
     P2=u[3];
           
     cosQ2= COS(Q2);
     cos2Q2= COS(2*Q2);

     aux1=-SIN(Q1);
     aux2=-SIN(Q2);
     aux3=2*SIN(2*Q2);
     aux4=-SIN(Q1+Q2);
     aux7=(C4-C5*cos2Q2);
     aux5=(C3*P1*P2)/aux7;
     aux6=(C5*(C2*P1*P1-C3*cosQ2*P1*P2+C1*P2*P2))/(aux7*aux7);

     f[0]= -((2*C2*P1-C3*cosQ2*P2)/aux7);
     f[1]= -((-C3*cosQ2*P1+2*C1*P2)/aux7);
     f[2]=aux1*C6+aux4*C7;
     f[3]=-aux3*aux6+aux4*C7-aux5*aux2-k*Q2;    


//     printf("\nOdefun4\n");
//     printf("%lg,%lg,%lg,%lg\n", f[0],f[1],f[2],f[3]);           

                        
     return;

}

val_type Ham4 (const int neq, const solution *u, const parameters *params)
{

/* ---------- First initializations -----------------------------------------*/

     val_type *uu;
     uu=u->uu;

/*------ declarations -------------------------------------------------*/ 

     val_type C1,C2,C3,C4,C5,C6,C7,k;
     val_type Q1,Q2,P1,P2;
     val_type cosQ1,cosQ2,cos2Q2,cosQ2Q1;
     val_type H;


/* ----------- implementation  ---------------------------------------*/


     C1=params->rpar[0];
     C2=params->rpar[1];
     C3=params->rpar[2];
     C4=params->rpar[3];
     C5=params->rpar[4];
     C6=params->rpar[5];
     C7=params->rpar[6];
     k=params->rpar[7];

     Q1=uu[0];
     Q2=uu[1];
     P1=uu[2];
     P2=uu[3];

     cosQ1= COS(Q1);
     cosQ2= COS(Q2);
     cos2Q2=COS(2*Q2);
     cosQ2Q1=COS(Q2+Q1);  

     H= (k*Q2*Q2/2-(C1*P2*P2+C2*P1*P1-C3*P2*P1*cosQ2) / (C4-C5*cos2Q2))
        - C6*cosQ1-C7*cosQ2Q1;

//     printf("\nHamiltonian\n");
//     printf("H=%lg\n",H);
     return(H);

}

/*------------------------------------------------------------------------------*/
/*										*/
/*        Problem-5:   Double PEndulum NEW 13.12.2016				*/
/*	    Ode5()=								*/
/*	    Ham5()=								*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode5 (const int neq,const val_type t,const val_type *u,
           val_type *f, const parameters *params)
{

 /* ---------- First initializations -----------------------------------------*/

 /*------ declarations -------------------------------------------------*/ 

     val_type C1,C2,C3,C4,C5,C6,C7;
     val_type Q1,Q2,P1,P2;
     val_type P2P1, cosQ2, cos2Q2;
     val_type sinQ1,sinQ2,sin2Q2,sinQ1Q2;
     val_type aux1,aux2;
     val_type daux1,daux2,daux3,daux4;

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
           

     P2P1=P2-P1;
     cosQ2= COS(Q2);
     cos2Q2= COS(2*Q2);
     sinQ1= SIN(Q1);
     sinQ2= SIN(Q2);
     sin2Q2= SIN(2*Q2);
     sinQ1Q2= SIN(Q1+Q2);
     
     aux1=C4-C5*cos2Q2;
     aux2=C3*cosQ2;

     daux1=(aux2*P2+2*C2*P2P1)/aux1;
     daux2=(2*C1*P2+aux2*P2P1)/aux1;
     daux3=(C3*P2*P2P1)/aux1;
     daux4=C5*(C1*P2*P2+aux2*P2*P2P1+C2*P2P1*P2P1)/(aux1*aux1);
    

     f[0]= -daux1;
     f[1]= daux1+daux2;    
     f[2]= -C7*sinQ1Q2-C6*sinQ1;
     f[3]= sinQ2*daux3+2*sin2Q2*daux4-C7*sinQ1Q2;


/*
     printf("\nOdefun5\n");
     printf("%.20lg,%.20lg,%.20lg,%.20lg\n", f[0],f[1],f[2],f[3]);           
*/
                        
     return;

}

val_type Ham5 (const int neq, const solution *u, const parameters *params)
{

/* ---------- First initializations -----------------------------------------*/

     val_type *uu;
     uu=u->uu;

/*------ declarations -------------------------------------------------*/ 

     val_type C1,C2,C3,C4,C5,C6,C7;
     val_type Q1,Q2,P1,P2;
     val_type P2P1, cosQ1,cosQ2, cos2Q2, cosQ1Q2;
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

     P2P1=P2-P1;
     cosQ1= COS(Q1);
     cosQ2= COS(Q2);
     cos2Q2= COS(2*Q2);
     cosQ1Q2= COS(Q1+Q2);

     H=((C1*P2*P2+C2*P2P1*P2P1+C3*P2*P2P1*cosQ2) / (C4-C5*cos2Q2))
        - C6*cosQ1-C7*cosQ1Q2;

//     printf("\nHamiltonian\n");
//     printf("H=%lg\n",H);
     return(H);

}


/*------------------------------------------------------------------------------*/
/*										*/
/*        Problem-6:	 							*/
/*	    Ode6()=								*/
/*	    Ham6()=								*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode6 (const int neq,const val_type t,const val_type *u,
           val_type *f, const parameters *params)

{
     return;
}

val_type Ham6 (const int neq, const solution *u, const parameters *params)
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

void Ode7 (const int neq,const val_type t,const val_type *u,
           val_type *f, const parameters *params)
{
     return;
}

val_type Ham7 (const int neq, const solution *u, const parameters *params)
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

void Ode8 (const int neq,const val_type t,const val_type *u,
           val_type *f, const parameters *params)
{
     return;
}

val_type Ham8 (const int neq, const solution *u, const parameters *params)
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

void Ode9 (const int neq,const val_type t,const val_type *u,
           val_type *f, const parameters *params)
{
     return;
}

val_type Ham9 (const int neq, const solution *u, const parameters *params)
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

void Ode10 (const int neq,const val_type t,const val_type *u,
           val_type *f, const parameters *params)
{
     return;
}

val_type Ham10 (const int neq, const solution *u, const parameters *params)
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


void Ode11 (const int neq,const val_type t,const val_type *u,
           val_type *f, const parameters *params)
{
 /* --------------------------------------------------------------------------*/
 /* specific Odependulum() for "ideal integrator" esperiment                  */
 /*     - input:      quadruple---> double                                    */
 /*     - operations: double                                                  */
 /*     - output:     double ---> quadruple.                                  */
 /* --------------------------------------------------------------------------*/

 /* ---------- First initializations -----------------------------------------*/

 /*------ declarations -------------------------------------------------------*/ 

     int i;

     double fdouble[neq];

     double C1,C2,C3,C4,C5,C6,C7;
     double Q1,Q2,P1,P2,Q12;
     double cosQ1Q2,sinQ1Q2,sinQ1,sinQ2;
     double aux0,aux1,dQ1aux1,aux2,aux3,aux4;

 /* ----------- implementation  ----------------------------------------------*/

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

     Q12=(Q1-Q2);
            
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
     fdouble[3]=(aux4-C7*sinQ2);            
              
     for (i=0; i<neq; i++) f[i]=fdouble[i];
          
     return ;

}



/*------------------------------------------------------------------------------*/
/*										*/
/*       NBody Problem:								*/
/*	       Ode12=OdeNBodyX (Ideal integrator)				*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode12 (const int neq,const val_type t,const val_type *u,
           val_type *f, const parameters *params)

{

 /* --------------------------------------------------------------------------*/
 /* specific Odependulum() for "ideal integrator" esperiment                  */
 /*     - input:      quadruple---> double                                    */
 /*     - operations: double                                                  */
 /*     - output:     double ---> quadruple.                                  */
 /* --------------------------------------------------------------------------*/


/* ---------- First initializations ------------------------------------------*/
  
     int dim,nd,nbody;
     dim=3;
     nbody=neq/(2*dim);
     nd=neq/2;
    
/*------ declarations --------------------------------------------------------*/ 

     int i,id,i1,i2,j,j1,j2;
     int ix,jx;
     double d3,qij;
     val_type *Gm;

     double fdouble[neq],Udouble[neq],Gmdouble[nbody];

/* ----------- implementation  -----------------------------------------------*/

/*   params->ipar=orderplanet (ascending order)*/

     Gm=params->rpar;

     for (i=0;i<nbody;i++) Gmdouble[i]=Gm[i];
     for (i=0;i<neq;i++) Udouble[i]=u[i];
    
     switch (params->eval)
     {
     case 1: /*OdePlanetsq*/

           for (ix=0; ix<nbody; ix++)
           {
                 i=params->ipar[ix];
                 i1=i*dim;
                 i2=nd+i1;
                 for (id=0; id<dim; id++) fdouble[i1+id]=Udouble[i2+id];
           }

           for (i=0; i<neq/2; i++) f[i]=fdouble[i];
              
     break;

     case 2: /*OdePlanetsv*/      

           for (ix=0; ix<nbody; ix++)
           {
                 i=params->ipar[ix];
                 i1=i*dim;
                 i2=nd+i1;
                 for (id=0; id<dim; id++) fdouble[i2+id]=0.;
           }
        
           for (ix=0; ix<nbody; ix++)
           {
                 i=params->ipar[ix];
                 i1=i*dim;
                 i2=nd+i1;

                 for (jx=ix+1; jx<nbody; jx++)
                 {
                       j=params->ipar[jx];
                       j1=j*dim;
                       j2=nd+j1;
                       d3=0.;
                       for (id=0; id<dim; id++) d3+=(u[i1+id]-u[j1+id])*(u[i1+id]-u[j1+id]);
                       d3=sqrt(d3)*d3;
//                       d3=pow(d3,3./2);

                       for (id=0; id<dim; id++) 
                       {
                             qij=(Udouble[i1+id]-Udouble[j1+id])/d3;
                             fdouble[i2+id]-=Gmdouble[j]*qij;
                             fdouble[j2+id]+=Gmdouble[i]*qij;
                       }   
                 }
           }
    
           for (i=neq/2; i<neq; i++) f[i]=fdouble[i];
     
     break;
                                     
     default:     
        
           for (ix=0; ix<nbody; ix++)
           {
                 i=params->ipar[ix];
                 i1=i*dim;
                 i2=nd+i1;
                 for (id=0; id<dim; id++) 
                 {
                       fdouble[i1+id]=Udouble[i2+id];
                       fdouble[i2+id]=0.;
                 }
           }

           for (ix=0; ix<nbody; ix++)
           {
		 i=params->ipar[ix];
                 i1=i*dim;
                 i2=nd+i1;
                 for (jx=ix+1; jx<nbody; jx++)
                 {
                       j=params->ipar[jx];
                       j1=j*dim;
                       j2=nd+j1;
                       d3=0.;
                       for (id=0; id<dim; id++) d3+=(u[i1+id]-u[j1+id])*(u[i1+id]-u[j1+id]);
                       d3=sqrt(d3)*d3;
//		       d3=pow(d3,3./2);

                       for (id=0; id<dim; id++) 
                       {
                             qij=(Udouble[i1+id]-Udouble[j1+id])/d3;
                             fdouble[i2+id]-=Gmdouble[j]*qij;
                             fdouble[j2+id]+=Gmdouble[i]*qij;
                       }   
                 }

           }

           for (i=0; i<neq; i++) f[i]=fdouble[i];
 
     break;

    }

    

    return ;

}

/*------------------------------------------------------------------------------*/
/*										*/
/*       OdePendulumStiff Problem: 						*/
/*		  Ode14= OdePendulumStiffX (Ideal integrator)			*/
/*										*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode14 (const int neq,const val_type t,const val_type *u,
           val_type *f, const parameters *params)

{
 /* --------------------------------------------------------------------------*/
 /* specific OdependulumStiff() for "ideal integrator" esperiment             */
 /*     - input:      quadruple---> double                                    */
 /*     - operations: double                                                  */
 /*     - output:     double ---> quadruple.                                  */
 /* --------------------------------------------------------------------------*/

 /* ---------- First initializations -----------------------------------------*/

 /*------ declarations -------------------------------------------------------*/ 

     int i;

     double fdouble[neq];

     double C1,C2,C3,C4,C5,C6,C7,k;
     double Q1,Q2,P1,P2;
     double cosQ2,cos2Q2;
     double aux1,aux2,aux3,aux4,aux5,aux6,aux7;

 /* ----------- implementation  ----------------------------------------------*/

     C1=params->rpar[0];
     C2=params->rpar[1];
     C3=params->rpar[2];
     C4=params->rpar[3];
     C5=params->rpar[4];
     C6=params->rpar[5];
     C7=params->rpar[6];
     k=params->rpar[7];

     Q1=u[0];
     Q2=u[1];
     P1=u[2];
     P2=u[3];
           
     cosQ2= cos(Q2);
     cos2Q2= cos(2*Q2);

     aux1=-sin(Q1);
     aux2=-sin(Q2);
     aux3=2*sin(2*Q2);
     aux4=-sin(Q1+Q2);
     aux7=(C4-C5*cos2Q2);
     aux5=(C3*P1*P2)/aux7;
     aux6=(C5*(C2*P1*P1-C3*cosQ2*P1*P2+C1*P2*P2))/(aux7*aux7);

     fdouble[0]= -((2*C2*P1-C3*cosQ2*P2)/aux7);
     fdouble[1]= -((-C3*cosQ2*P1+2*C1*P2)/aux7);
     fdouble[2]=aux1*C6+aux4*C7;
     fdouble[3]=-aux3*aux6+aux4*C7-aux5*aux2-k*Q2;   

     for (i=0; i<neq; i++) f[i]=fdouble[i];

    return ;

}


/*------------------------------------------------------------------------------*/
/*										*/
/*       OdePendulum Problem: 	     NEW 13.12.2016				*/
/*	 Ode15= OdePendulumX (Ideal integrator)	                 		*/
/*										*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode15 (const int neq,const val_type t,const val_type *u,
           val_type *f, const parameters *params)

{
 /* --------------------------------------------------------------------------*/
 /* specific Odependulum() for "ideal integrator" esperiment                  */
 /*     - input:      quadruple---> double                                    */
 /*     - operations: double                                                  */
 /*     - output:     double ---> quadruple.                                  */
 /* --------------------------------------------------------------------------*/

 /* ---------- First initializations -----------------------------------------*/

 /*------ declarations -------------------------------------------------------*/ 

     int i;

     double fdouble[neq];

     double C1,C2,C3,C4,C5,C6,C7;
     double Q1,Q2,P1,P2;
     double P2P1, cosQ2, cos2Q2;
     double sinQ1,sinQ2,sin2Q2,sinQ1Q2;
     double aux1,aux2;
     double daux1,daux2,daux3,daux4;

 /* ----------- implementation  ----------------------------------------------*/

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
           
     P2P1=P2-P1;
     cosQ2= COS(Q2);
     cos2Q2= COS(2*Q2);
     sinQ1= SIN(Q1);
     sinQ2= SIN(Q2);
     sin2Q2= SIN(2*Q2);
     sinQ1Q2= SIN(Q1+Q2);
     
     aux1=C4-C5*cos2Q2;
     aux2=C3*cosQ2;

     daux1=(aux2*P2+2*C2*P2P1)/aux1;
     daux2=(2*C1*P2+aux2*P2P1)/aux1;
     daux3=(C3*P2*P2P1)/aux1;
     daux4=C5*(C1*P2*P2+aux2*P2*P2P1+C2*P2P1*P2P1)/(aux1*aux1);
    
     fdouble[0]= -daux1;
     fdouble[1]= daux1+daux2;    
     fdouble[2]= -C7*sinQ1Q2-C6*sinQ1;
     fdouble[3]= sinQ2*daux3+2*sin2Q2*daux4-C7*sinQ1Q2;

     for (i=0; i<neq; i++) f[i]=fdouble[i];

    return ;

}




