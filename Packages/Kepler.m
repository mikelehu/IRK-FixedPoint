(* ::Package:: *)

(* ::Section:: *)
(*Kepler Problem *)


BeginPackage["Kepler`"];

KeplerODE::usage="arameters[[-1]]=1 ----> Return (fq)
    parameters[[-1]]=2 ----> Return (fp)
    otherwise        ----> Return (f=fq,fp) ";
KeplerODEc::usage=" ";
KeplerODEd::usage=" ";
KeplerODEfunc::usage=" ";
KeplerHam::usage=" 
    H(q,v)= ";

Begin["`Private`"];


(* ::Subsection:: *)
(*Functions*)


(* ::Subsubsection::Closed:: *)
(*Odefun*)


KeplerODE[t_,u_,parameters_List]:=
Module[{q1,q2,p1,p2,r3},
q1=u[[1]];
q2=u[[2]];
p1=u[[3]];
p2=u[[4]];

If[parameters[[-1]]==1,
{p1,p2},
r3=(q1^2+q2^2)^(3/2);
If[parameters[[-1]]==2,
{-q1/r3,-q2/r3},
{p1,p2,-q1/r3,-q2/r3}]
]
];

KeplerODEc[t_,u_,parameters_List]:=Apply[KeplerODEfunc,Flatten[{t,u,parameters}]];
KeplerODEd[t_,u_,parameters_List]:=SetPrecision[Apply[KeplerODEfunc,Flatten[N[{t,u,parameters}]]],40];

KeplerODEfunc=Compile[{t,q1,q2,p1,p2,par1},
Module[{r3},
If[par1==1,
{p1,p2},
r3=(q1^2+q2^2)^(3/2);
If[par1==2,
{-q1/r3,-q2/r3},
{p1,p2,-q1/r3,-q2/r3}]
]

]
];



(* ::Input:: *)
(**)


(* ::Subsubsection::Closed:: *)
(*Hamiltonian*)


KeplerHam[{t_,u_,e_},parameters_List, doi_:100]:=
Module[{q1,q2,p1,p2},
uu=SetPrecision[u,doi]+SetPrecision[e,doi];
q1=uu[[1]];
q2=uu[[2]];
p1=uu[[3]];
p2=uu[[4]];
(1/2)*(p1^2+p2^2)-1/Sqrt[q1^2+q2^2]
];



(* ::Subsection:: *)
(*End*)


End[];
EndPackage[];

