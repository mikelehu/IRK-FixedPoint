/*----------------------------------------------------------------------------*/
/*									      */
/*         GaussInitData.c						      */
/*									      */
/*         Functions: InitialDta()					      */
/*         Initial values:						      */
/*		1:  Double-Pendulum non-Chaotic.			      */
/*		2:  Double-Pendulum Chaotic.				      */
/*		3: Outer-Solar-System).				      */
/* ---------------------------------------------------------------------------*/


#include <GaussInitData.h>


void InitialData (toptions *options, solution *u,
                  parameters *params,ode_sys *system)

{

     int i;
     val_type m1,m2,l1,l2,g;
     val_type C1,C2,C3,C4,C5,C6,C7;
     val_type GK;
    
     switch (system->problem)
     {
     case 1:     /* Double Pendulum: non-Chaotic initial values. */

           system->neq=4;
           system->n=2;
           options->h = POW(2,-7);
           options->t0=0.;
           options->t1=pow(2,12);  
           options->sampling=POW(2,10);

           m1=1.;
           m2=1.;
           l1=1.;
           l2=1.;
           g=9.8;

           //u0=(q,p) 
           u->uu[0]=11./10; u->uu[1]=-11./10; 	   
           u->uu[2]=13873./5000;  u->uu[3]=13873./5000;

           //e0
           u->ee[0]=-8.8817841970012523234e-17;
           u->ee[1]= 8.8817841970012523234e-17;
           u->ee[2]=4.4764192352886311710e-17; 
           u->ee[3]=4.4764192352886311710e-17;      
       	
           C1=-l1*l1*(m1+m2);
           C2=-l2*l2*m2;
           C3=-2*l1*l2*m2;
           C4=-2*l1*l1*l2*l2*m2*m1-l1*l1*l2*l2*m2*m2;
           C5=-l1*l1*l2*l2*m2*m2;
           C6=g*l1*(m1+m2);
           C7=g*l2*m2;

           params->rpar[0]=C1;
           params->rpar[1]=C2;	
           params->rpar[2]=C3;
           params->rpar[3]=C4;
           params->rpar[4]=C5;
           params->rpar[5]=C6;
           params->rpar[6]=C7;	 

     break;	

     case 2:     /* Double Pendulum: Chaotic initial values. */

           system->neq=4;
           system->n=2;
           options->h = POW(2,-7);
           options->t0=0.;
           options->t1=pow(2,8);    
           options->sampling=POW(2,6);            	  

           m1=1.;
           m2=1.;
           l1=1.;
           l2=1.;
           g=9.8;

           //u0=(q,p) 
           u->uu[0]=0.; u->uu[1]=0.; 	   
           u->uu[2]=3873./1000; u->uu[3]=3873./1000;

           //e0
           u->ee[0]=0.;
           u->ee[1]=0.;
           u->ee[2]=-2.2026824808563105762e-16;
           u->ee[3]=-2.2026824808563105762e-16;

           C1=-l1*l1*(m1+m2);
           C2=-l2*l2*m2;
           C3=-2*l1*l2*m2;
           C4=-2*l1*l1*l2*l2*m2*m1-l1*l1*l2*l2*m2*m2;
           C5=-l1*l1*l2*l2*m2*m2;
           C6=g*l1*(m1+m2);
           C7=g*l2*m2;

           params->rpar[0]=C1;
           params->rpar[1]=C2;	
           params->rpar[2]=C3;
           params->rpar[3]=C4;
           params->rpar[4]=C5;
           params->rpar[5]=C6;
           params->rpar[6]=C7;

     break;	

     case 3:      /*Outer-Solar-System)*/

           system->neq=36;
           system->n=6;
           options->h = 500./3;
           options->t0=0;
           options->t1=10000000.;
           options->sampling=120;

           // Gm
           GK=2.95912208286*pow(10,-4);
           params->rpar[0]=GK*1.00000597682;		// sun+inner planets
           params->rpar[1]=GK*0.000954786104043;	// Jupiter
           params->rpar[2]=GK*0.000285583733151;	// Saturn
           params->rpar[3]=GK*0.0000437273164546;	// Uranus
           params->rpar[4]=GK*0.0000517759138449;	// Neptune
           params->rpar[5]=GK*1.0/(1.3*pow(10,8));	// Pluto

           // orderplanets
           params->ipar[0]=5;
           params->ipar[1]=3;
	   params->ipar[2]=4;
	   params->ipar[3]=2;
           params->ipar[4]=1;
           params->ipar[5]=0;

           //position
 
    	   u->uu[0]=-2.047098298789102e-4; u->uu[1]=6.550139855052494e-3;	u->uu[2]=2.824833990245128e-3;
           u->uu[3]=-3.502570009829879;	   u->uu[4]=-3.810434560144947;		u->uu[5]=-1.547971466009755;
           u->uu[6]= 9.07532669017012;	   u->uu[7]=-3.039285160144947;		u->uu[8]=-1.645545966009755;
           u->uu[9]= 8.30993729017012;	   u->uu[10]=-1.628355846014495e1;	u->uu[11]=-7.249302966009755;
           u->uu[12]=1.147056189017012e1;  u->uu[13]=-2.572293276014495e1;	u->uu[14]=-1.081412076600976e1;
           u->uu[15]=-1.553894040982988e1; u->uu[16]=-2.521600926014495e1;	u->uu[17]=-3.187413366009755;
      
           //velocities

           u->uu[18]=-6.175529636225841e-6;	 u->uu[19]=2.43502570182194e-6;	        u->uu[20]=1.223839570932369e-6;
           u->uu[21]=5.648114470363774e-3;       u->uu[22]=-4.122464974298178e-3;	u->uu[23]=-1.904666160429068e-3;
           u->uu[24]=1.677004470363774e-3;       u->uu[25]=4.837685025701822e-3;	u->uu[26]=1.925843839570932e-3;
           u->uu[27]=3.535604470363774e-3;       u->uu[28]=1.373455025701822e-3;	u->uu[29]=5.515138395709323e-4;
           u->uu[30]=2.883124470363774e-3;       u->uu[31]=1.147705025701822e-3;	u->uu[32]=3.979938395709323e-4;
           u->uu[33]=2.761074470363774e-3;       u->uu[34]=-1.704584974298178e-3;	u->uu[35]=-1.363816160429068e-3;

           for (i=0; i<system->neq; i++) { u->ee[i]=0.;}

     break;


     default:
           printf("Initialdata.c: problem=%i does not exist\n",
                   system->problem);
     break;
     }



}


