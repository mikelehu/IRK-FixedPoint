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
/*										*/
/*------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <def.h>


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
     printf("\nOdefun1\n");
     printf("%.20lg,%.20lg,%.20lg,%.20lg\n", f[0],f[1],f[2],f[3]);           
*/
                        
     return;
}



val_type Ham1 (const int neq, const solution *u, const parameters *params)
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

/*
     printf("\nHamiltonian\n");
     printf("H=%lg\n",H);
*/
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
                        
     return;

}

val_type Ham4 (const int neq, const solution *u, const parameters *params)
{

     return(0.);

}

/*------------------------------------------------------------------------------*/
/*										*/
/*        Problem-5:   								*/
/*	    Ode5()=								*/
/*	    Ham5()=								*/
/*										*/
/*------------------------------------------------------------------------------*/

void Ode5 (const int neq,const val_type t,const val_type *u,
           val_type *f, const parameters *params)
{

                        
     return;

}

val_type Ham5 (const int neq, const solution *u, const parameters *params)
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


