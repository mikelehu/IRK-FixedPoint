(* ::Package:: *)

(* ::Section:: *)
(*Double Pendulum  Problem (NEW 13-12-2016) *)


BeginPackage["DoublePendulumNEW`"];

DoublePendulumODENEW::usage=" ";
DoublePendulumHamNEW::usage="H(q,v)= ";
DoublePendulumJACNEW::usage="";    

DoublePendulumHam3NEW::usage="Fortran ";


Begin["`Private`"];


(* ::Subsection:: *)
(*Functions*)


(* ::Subsubsection::Closed:: *)
(*Odefun*)


DoublePendulumODENEW[t_,u_,parameters_List]:=
Module[{\[Phi], \[Theta],p\[Phi],p\[Theta],
C1,C2,C3,C4,C5,C6,C7,
p\[Theta]p\[Phi],cos\[Theta],cos2\[Theta],
sin\[Phi],sin\[Theta],sin2\[Theta],sin\[Theta]\[Phi],
aux1,aux2,
daux1,daux2,daux3,daux4
},

\[Phi]=u[[1]];
\[Theta]=u[[2]];
p\[Phi]=u[[3]];
p\[Theta]=u[[4]];

C1=parameters[[1]];
C2=parameters[[2]];
C3=parameters[[3]];
C4=parameters[[4]];
C5=parameters[[5]];
C6=parameters[[6]];
C7=parameters[[7]];


p\[Theta]p\[Phi]=p\[Theta]-p\[Phi];
cos\[Theta]=Cos[\[Theta]];
cos2\[Theta]=Cos[2 \[Theta]];
sin\[Phi]=Sin[\[Phi]];
sin\[Theta]=Sin[\[Theta]];
sin2\[Theta]=Sin[2 \[Theta]];
sin\[Theta]\[Phi]=Sin[\[Theta]+\[Phi]];

aux1=C4-C5 cos2\[Theta];
aux2=C3 cos\[Theta] ;

daux1=(aux2 p\[Theta]+2 C2 p\[Theta]p\[Phi])/aux1;
daux2=(2 C1 p\[Theta]+aux2 p\[Theta]p\[Phi])/aux1;
daux3=(C3 p\[Theta] p\[Theta]p\[Phi])/aux1;
daux4=(C5 (C1 p\[Theta]^2+aux2 p\[Theta] p\[Theta]p\[Phi]+C2 p\[Theta]p\[Phi]^2))/(aux1*aux1);

{-daux1,
daux1+daux2,
-C7 sin\[Theta]\[Phi] -C6 sin\[Phi],
sin\[Theta] daux3+2 sin2\[Theta] daux4-C7 sin\[Theta]\[Phi]
}
];



(* ::Subsubsection::Closed:: *)
(*Hamiltonian-C*)


DoublePendulumHamNEW[neq_,out_,parameters_List, doi_:100]:=
Module[{uu,
 \[Phi], \[Theta],p\[Phi],p\[Theta],
c1,c2,c3,c4,c5,c6,c7,
p\[Theta]p\[Phi],cos\[Phi],cos\[Theta],cos2\[Theta],cos\[Theta]\[Phi]},

uu=SetPrecision[out[[2;;5]],doi]+SetPrecision[out[[6;;9]],doi];
\[Phi]=uu[[1]];
\[Theta]=uu[[2]];
p\[Phi]=uu[[3]];
p\[Theta]=uu[[4]];

c1=parameters[[1]]; 
c2=parameters[[2]]; 
c3=parameters[[3]]; 
c4=parameters[[4]]; 
c5=parameters[[5]]; 
c6=parameters[[6]]; 
c7=parameters[[7]];  

p\[Theta]p\[Phi]=(p\[Theta]-p\[Phi]);
cos\[Phi]=Cos[\[Phi]];
cos\[Theta]= Cos[ \[Theta]];
cos2\[Theta]=Cos[2  \[Theta]];
cos\[Theta]\[Phi]=Cos[\[Theta]+\[Phi]];

(c1 p\[Theta]^2+c2 p\[Theta]p\[Phi]^2+c3  p\[Theta] p\[Theta]p\[Phi] cos\[Theta])/(c4-c5 cos2\[Theta])-c6 cos\[Phi]-c7 cos\[Theta]\[Phi]

];


(* ::Subsubsection:: *)
(*Jacobian*)


DoublePendulumJACNEW[t_,u_,parameters_List]:=
Module[{c1,c2,c3,c4,c5,c6,c7,
         Q1,Q2,P1,P2,    
P2P1,cosQ1,cosQ2,cos2Q2,cosQ1Q2,sinQ1,sinQ2,sin2Q2,sinQ1Q2,
aux1,aux2,aux3,aux4,aux5,aux6,
aux11,aux12,aux13,aux21,aux22,aux23,aux24,aux41,aux42,aux43,aux44,
aux51,aux52,aux53,
daux1,daux41,daux42,
aux101,aux102,aux103,aux104,aux105,aux106,aux107,aux108,aux109,
aux110,aux111,aux112,aux113,aux114},

c1=parameters[[1]];
c2=parameters[[2]];
c3=parameters[[3]];
c4=parameters[[4]];
c5=parameters[[5]];
c6=parameters[[6]];
c7=parameters[[7]];

Q1=u[[1]];
       Q2=u[[2]];
       P1=u[[3]];
       P2=u[[4]];

P2P1=-P1+P2;
cosQ1=Cos[Q1];
cosQ2=Cos[Q2];
cos2Q2=Cos[2 Q2];
cosQ1Q2=Cos[Q1+Q2];
sinQ1=Sin[Q1];
sinQ2=Sin[Q2];
sin2Q2=Sin[2 Q2];
sinQ1Q2=Sin[Q1+Q2];

aux1=c4-c5 cos2Q2;
aux2=c3 cosQ2;
aux3=c6 cosQ1;
aux4=c7 cosQ1Q2;
aux5=2c5 sin2Q2;
aux6=-c3 sinQ2;

aux11=aux2 P2;
aux12=aux2 P2P1;
aux13=aux2 P2 P2P1;

aux21=aux6 P2;
aux22=aux6 P2P1;
aux23=aux6 P2 P2P1;
aux24=-2aux23 aux5 ;

aux41=aux1*aux1;
aux42=aux1*aux41;
aux43=P2*P2;
aux44=P2P1*P2P1;

aux51=c1 aux43+aux13 +c2 aux44;
aux52=( aux11+2 c2 P2P1);
aux53=2 c1 P2+aux12;

daux1=1/aux1;
daux41=1/aux41;
daux42=1/aux42;

aux101=daux1 aux21 -aux5 aux52 daux41;
aux102=-2 c2 daux1;
aux103=(aux2+2 c2)daux1;
aux104=aux22 daux1-aux5 aux53 daux41;
aux105=2(aux2+ c1+ c2) daux1;
aux106=-c3 P2 daux1 sinQ2;
aux107=(c3 P2+c3 P2P1) sinQ2 daux1;
aux108=2 c5 aux5 aux51 daux42;
aux109= aux13 daux1;
aux110=4 c5 aux51 cos2Q2 daux41;
aux111=aux24 daux41;
aux112=- 2 c5 aux52 sin2Q2 daux41;
aux113=2c5 daux41 (aux53+aux52)sin2Q2;
aux114=2sin2Q2 aux108;


{{-aux3 -aux4,-aux4,0,0},
{-aux4,
aux110+aux109-aux4-aux111-aux114,
aux106 +aux112 ,aux107 +aux113 },
{0,-aux101,-aux102,-aux103},
{0,aux101+aux104,-aux103,aux105}}
];





(* ::Subsubsection::Closed:: *)
(*Hamiltonian Fortran*)


DoublePendulumHam3NEW[neq_,out_,parameters_List, doi_:100]:=
Module[{uu,
 \[Phi], \[Theta],p\[Phi],p\[Theta],
c1,c2,c3,c4,c5,c6,c7,
p\[Theta]p\[Phi],cos\[Phi],cos\[Theta],cos2\[Theta],cos\[Theta]\[Phi]},

(*uu=SetPrecision[out[[2;;5]],doi]+SetPrecision[out[[6;;9]],doi];*)
uu=SetPrecision[out[[2;;5]],doi];
\[Phi]=uu[[1]];
\[Theta]=uu[[2]];
p\[Phi]=uu[[3]];
p\[Theta]=uu[[4]];

c1=parameters[[1]]; 
c2=parameters[[2]]; 
c3=parameters[[3]]; 
c4=parameters[[4]]; 
c5=parameters[[5]]; 
c6=parameters[[6]]; 
c7=parameters[[7]];  

p\[Theta]p\[Phi]=p\[Theta]-p\[Phi];
cos\[Phi]=Cos[\[Phi]];
cos\[Theta]= Cos[ \[Theta]];
cos2\[Theta]=Cos[2  \[Theta]];
cos\[Theta]\[Phi]=Cos[\[Theta]+\[Phi]];

(c1 p\[Theta]^2+c2 p\[Theta]p\[Phi]^2+c3  p\[Theta] p\[Theta]p\[Phi] cos\[Theta])/(c4-c5 cos2\[Theta])-c6 cos\[Phi]-c7 cos\[Theta]\[Phi]

];


(* ::Subsection:: *)
(*End*)


End[];
EndPackage[];



