(* ::Package:: *)

(* ::Section:: *)
(*Double Pendulum Problem *)


BeginPackage["DoublePendulum`"];

DoublePendulumODE::usage=" ";
DoublePendulumODEc::usage=" ";
DoublePendulumODEd::usage=" ";
DoublePendulumODEfunc::usage=" ";

DoublePendulumHam2::usage="H(q,v)= ";
DoublePendulumJAC::usage="";
DoublePendulumJACc::usage=" ";
DoublePendulumJACfunc::usage=" ";

Begin["`Private`"];


(* ::Subsection:: *)
(*Functions*)


(* ::Subsubsection::Closed:: *)
(*Odefun*)


DoublePendulumODE[t_,u_,parameters_List]:=
Module[{C1,C2,C3,C4,C5,C6,C7,
Q1,Q2,P1,P2,
Q12,cosQ1Q2,sinQ1Q2,sin2Q1Q2,sinQ1,sinQ2,cosQ1,cosQ2,
aux0,aux1,dQ1aux1,aux2,aux3,aux4},

Q1=u[[1]];
Q2=u[[2]];
P1=u[[3]];
P2=u[[4]];

C1=parameters[[1]];
C2=parameters[[2]];
C3=parameters[[3]];
C4=parameters[[4]];
C5=parameters[[5]];
C6=parameters[[6]];
C7=parameters[[7]];

Q12 = Q1-Q2;
cosQ1Q2=Cos[Q12];
sinQ1Q2=Sin[Q12];
sin2Q1Q2=Sin[2Q12];
sinQ1=Sin[Q1];
sinQ2=Sin[Q2];
cosQ1=Cos[Q1];
cosQ2=Cos[Q2];

 aux0=C5 (sinQ1Q2)^2;
aux1= C4+aux0;
dQ1aux1=2*C5*sinQ1Q2*cosQ1Q2;
aux2=(aux1)^2;
aux3= C3*cosQ1Q2;
aux4=(-1/aux2)*(C1*P1^2+C2*P2^2+P1*P2*aux3)*dQ1aux1-(C3*P1*P2*sinQ1Q2)/aux1;

{(2*C1*P1+aux3*P2)/aux1,
(2*C2*P2+aux3*P1)/aux1,
-(aux4+C6*sinQ1),
-(-aux4+C7*sinQ2)}];

DoublePendulumODEc[t_,u_,parameters_List]:=Apply[DoublePendulumODEfunc,Flatten[{t,u,parameters}]];
Quadprec=Abs[Log[10.,2^-113]];
DoublePendulumODEd[t_,u_,parameters_List]:=SetPrecision[Apply[DoublePendulumODEfunc,Flatten[N[{t,u,parameters}]]],Quadprec];

DoublePendulumODEfunc=Compile[{t,Q1,Q2,P1,P2,C1,C2,C3,C4,C5,C6,C7},
Module[{Q12,cosQ1Q2,sinQ1Q2,sin2Q1Q2,sinQ1,sinQ2,cosQ1,cosQ2,
aux0,aux1,dQ1aux1,aux2,aux3,aux4},

Q12 = Q1-Q2;
cosQ1Q2=Cos[Q12];
sinQ1Q2=Sin[Q12];
sin2Q1Q2=Sin[2Q12];
sinQ1=Sin[Q1];
sinQ2=Sin[Q2];
cosQ1=Cos[Q1];
cosQ2=Cos[Q2];

 aux0=C5 (sinQ1Q2)^2;
aux1= C4+aux0;
dQ1aux1=2*C5*sinQ1Q2*cosQ1Q2;
aux2=(aux1)^2;
aux3= C3*cosQ1Q2;
aux4=(-1/aux2)*(C1*P1^2+C2*P2^2+P1*P2*aux3)*dQ1aux1-(C3*P1*P2*sinQ1Q2)/aux1;

{(2*C1*P1+aux3*P2)/aux1,
(2*C2*P2+aux3*P1)/aux1,
-(aux4+C6*sinQ1),
-(-aux4+C7*sinQ2)}]];



(* ::Subsubsection::Closed:: *)
(*Hamiltonian C*)


DoublePendulumHam2[neq_,out_,parameters_List, doi_:100]:=
Module[{uu,C1,C2,C3,C4,C5,C6,C7,
Q1,Q2,P1,P2,
Q12,cosQ1Q2,sinQ1Q2,cosQ1,cosQ2},


uu=SetPrecision[out[[2;;5]],doi]+SetPrecision[out[[6;;9]],doi];
Q1=uu[[1]];
Q2=uu[[2]];
P1=uu[[3]];
P2=uu[[4]];

C1=parameters[[1]];
C2=parameters[[2]];
C3=parameters[[3]];
C4=parameters[[4]];
C5=parameters[[5]];
C6=parameters[[6]];
C7=parameters[[7]];


Q12=Q1-Q2;
cosQ1Q2=Cos[Q1-Q2];
sinQ1Q2=Sin[Q1-Q2];
cosQ1=Cos[Q1];
cosQ2=Cos[Q2];


((C1*P1^2+C2*P2^2+C3*P1*P2*cosQ1Q2)/(C4+C5*sinQ1Q2^2))-C6*cosQ1-C7*cosQ2];


(* ::Subsubsection:: *)
(*Jacobian*)


DoublePendulumJAC[t_,u_,parameters_List]:=
Module[{C1,C2,C3,C4,C5,C6,C7,
         Q1,Q2,P1,P2,
        Q12,cosQ1Q2,sinQ1Q2,cosQ1,cosQ2,
        aux0,aux1,dQ1aux1,aux2,aux3,
        daux10,daux104,daux17,daux22,daux33,daux34,daux37,daux38,   
        daux41,daux42,daux45,daux46,daux49,daux5,daux50,daux53,daux54,daux55,daux56,daux6,daux66,daux69,daux80,daux83,daux89,daux9},

C1=parameters[[1]];
C2=parameters[[2]];
C3=parameters[[3]];
C4=parameters[[4]];
C5=parameters[[5]];
C6=parameters[[6]];
C7=parameters[[7]];

Q1=u[[1]];
Q2=u[[2]];
P1=u[[3]];
P2=u[[4]];
Q12=Q1-Q2;

cosQ1Q2=Cos[Q12];
sinQ1Q2=Sin[Q12];
cosQ1=Cos[Q1];
cosQ2=Cos[Q2];

aux0=C5 sinQ1Q2^2;
aux1=aux0+C4;
dQ1aux1=2 C5 cosQ1Q2 sinQ1Q2;
aux2=aux1^2;
aux3=C3 cosQ1Q2;

daux5=-sinQ1Q2;
daux6=sinQ1Q2;
daux9=cosQ1Q2;
daux10=-cosQ1Q2; 
daux17=cosQ1;
daux22=cosQ2;
daux33=2 C5 daux9 sinQ1Q2;
daux34=2 C5 daux10 sinQ1Q2;
daux37=daux33;
daux38=daux34;
daux41=2 C5 cosQ1Q2 daux9+2 C5 daux5 sinQ1Q2;
daux42=2 C5 cosQ1Q2 daux10+2 C5 daux6 sinQ1Q2;
daux45=2 aux1 daux37;
daux46=2 aux1 daux38;
daux49=C3 daux5;

daux50=C3 daux6;
daux53=-((C3 daux9 P1 P2)/aux1)-(daux49 dQ1aux1 P1 P2)/aux2-
         (daux41 (C1 P1^2+aux3 P1 P2+C2 P2^2))/aux2+
         (daux45 dQ1aux1 (C1 P1^2+aux3 P1 P2+C2 P2^2))/aux2^2+(C3 daux37 P1 P2 sinQ1Q2)/aux1^2;
daux54=-((C3 daux10 P1 P2)/aux1)-(daux50 dQ1aux1 P1 P2)/aux2-
         (daux42 (C1 P1^2+aux3 P1 P2+C2 P2^2))/aux2+
         (daux46 dQ1aux1 (C1 P1^2+aux3 P1 P2+C2 P2^2))/aux2^2+(C3 daux38 P1 P2 sinQ1Q2)/aux1^2;
daux55=-((dQ1aux1 (2 C1 P1+aux3 P2))/aux2)-(C3 P2 sinQ1Q2)/aux1;
daux56=-((dQ1aux1 (aux3 P1+2 C2 P2))/aux2)-(C3 P1 sinQ1Q2)/aux1;

daux66=-((2 C1 P1+aux3 P2)/aux1^2);
daux69=P2/aux1;
daux80=-((aux3 P1+2 C2 P2)/aux1^2);
daux83=P1/aux1;
daux89=-C6;
daux104=-C7;

{{daux37 daux66+daux49 daux69,daux38 daux66+daux50 daux69,(2 C1)/aux1,aux3/aux1},
{daux37 daux80+daux49 daux83,daux38 daux80+daux50 daux83,aux3/aux1,(2 C2)/aux1},
{-daux53+daux17 daux89,-daux54,-daux55,-daux56},
{daux53,daux104 daux22+daux54,daux55,daux56}}
];



DoublePendulumJACc[t_,u_,parameters_List]:=Apply[DoublePendulumJACfunc,Flatten[{t,u,parameters}]];

DoublePendulumJACfunc=Compile[{t,Q1,Q2,P1,P2,C1,C2,C3,C4,C5,C6,C7},
Module[{Q12,cosQ1Q2,sinQ1Q2,sin2Q1Q2,sinQ1,sinQ2,cosQ1,cosQ2,
        aux0,aux1,dQ1aux1,aux2,aux3,aux4,
        daux10,daux104,daux13,daux14,daux17,daux22,daux25,daux30,daux33,daux34,daux37,daux38, 
        daux41,daux42,daux45,daux46,daux49,daux5,daux50,daux53,daux54,daux55,daux56,daux6,
        daux66,daux69,daux80,daux83,daux89,daux9},

Q12=Q1-Q2;
cosQ1Q2=Cos[Q12];
sinQ1Q2=Sin[Q12];
sin2Q1Q2=Sin[2 Q12];
sinQ1=Sin[Q1];
sinQ2=Sin[Q2];
cosQ1=Cos[Q1];
cosQ2=Cos[Q2];
aux0=C5 sinQ1Q2^2;
aux1=aux0+C4;
dQ1aux1=2 C5 cosQ1Q2 sinQ1Q2;
aux2=aux1^2;
aux3=C3 cosQ1Q2;
aux4=-((dQ1aux1 (C1 P1^2+aux3 P1 P2+C2 P2^2))/aux2)-(C3 P1 P2 sinQ1Q2)/aux1;
daux5=-Sin[Q12];
daux6=Sin[Q12];
daux9=Cos[Q12];
daux10=-Cos[Q12];daux104=-C7;
daux13=2 Cos[2 Q12];daux14=-2 Cos[2 Q12];
daux17=Cos[Q1];daux22=Cos[Q2];
daux25=-Sin[Q1];daux30=-Sin[Q2];
daux33=2 C5 daux9 sinQ1Q2;
daux34=2 C5 daux10 sinQ1Q2;
daux37=daux33;daux38=daux34;
daux41=2 C5 cosQ1Q2 daux9+2 C5 daux5 sinQ1Q2;
daux42=2 C5 cosQ1Q2 daux10+2 C5 daux6 sinQ1Q2;
daux45=2 aux1 daux37;
daux46=2 aux1 daux38;
daux49=C3 daux5;
daux50=C3 daux6;
daux53=-((C3 daux9 P1 P2)/aux1)-(daux49 dQ1aux1 P1 P2)/aux2-(daux41 (C1 P1^2+aux3 P1 P2+C2 P2^2))/aux2+(daux45 dQ1aux1 (C1 P1^2+aux3 P1 P2+C2 P2^2))/aux2^2+(C3 daux37 P1 P2 sinQ1Q2)/aux1^2;
daux54=-((C3 daux10 P1 P2)/aux1)-(daux50 dQ1aux1 P1 P2)/aux2-(daux42 (C1 P1^2+aux3 P1 P2+C2 P2^2))/aux2+(daux46 dQ1aux1 (C1 P1^2+aux3 P1 P2+C2 P2^2))/aux2^2+(C3 daux38 P1 P2 sinQ1Q2)/aux1^2;
daux55=-((dQ1aux1 (2 C1 P1+aux3 P2))/aux2)-(C3 P2 sinQ1Q2)/aux1;
daux56=-((dQ1aux1 (aux3 P1+2 C2 P2))/aux2)-(C3 P1 sinQ1Q2)/aux1;
daux66=-((2 C1 P1+aux3 P2)/aux1^2);
daux69=P2/aux1;
daux80=-((aux3 P1+2 C2 P2)/aux1^2);
daux83=P1/aux1;
daux89=-C6;

{{daux37 daux66+daux49 daux69,daux38 daux66+daux50 daux69,(2 C1)/aux1,aux3/aux1},
{daux37 daux80+daux49 daux83,daux38 daux80+daux50 daux83,aux3/aux1,(2 C2)/aux1},
{-daux53+daux17 daux89,-daux54,-daux55,-daux56},
{daux53,daux104 daux22+daux54,daux55,daux56}}
]];





(* ::Subsection:: *)
(*End*)


End[];
EndPackage[];



