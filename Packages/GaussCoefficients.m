(* ::Package:: *)

(* ::Section:: *)
(*GaussCoefficients*)


BeginPackage["GaussCoefficients`"];

NN::usage=" ";
GaussCoefficients::usage=" ";
InterpolationCoeff::usage=" ";

Begin["`Private`"];


(* ::Subsection:: *)
(*Functions*)


(* ::Subsubsection::Closed:: *)
(*Functions*)


NN[y_,doi_]:=Module[{d,exp,naux,m},
(*exp=Last[RealDigits[y,2]];*)
{d,exp}=RealDigits[y,2,doi];
naux=-exp+doi;
m=2^naux;
Round[y*m]/m
];
SetAttributes[NN,Listable];


(* ::Text:: *)
(*Float:    32 bit -->   1 + 8 + 23  -->  eps=2^(-24)  = 6e-8.*)
(*Double: 64 bit -->  1 +11 + 52 -->   eps=2^(-53) = 1e-16.*)
(*Quadruple: 128 bit --> 1 + 15 + 112 --> eps=2^(-113) = 1e-36. *)


GaussCoefficients[s_,prec_]:=
Module[{highprec,ks,i,j, bits,
        a,b,c,m,aD,bD,cD,mD,nuD,aDD,bDD,cDD,mDD},

highprec=200;

ks=Flatten[NDSolve`ImplicitRungeKuttaGaussCoefficients[2*s,highprec],1];
a=Simplify[Table[ks[[i]],{i,1,s}]];
b=Table[ks[[i+s]],{i,1,s}];
c=Table[ks[[i+2 s]],{i,1,s}];
m=Table[a[[i,j]]/b[[j]],{i,s},{j,s}];

bits=Round[Log[2.,10] prec];

mD=Map[NN[#,bits] &,m];
Do[Do[mD[[i,j]]=NN[mD[[i,j]],bits],{j,i-1}],{i,s}];
Do[Do[mD[[j,i]]=1-mD[[i,j]],{j,i-1}],{i,s}];

bD=Map[NN[#,bits] &,b];

aDD= mD . DiagonalMatrix[bD];
cDD=Table[Sum[mD[[i,j]]*bD[[j]],{j,s}],{i,s}]; 

aD=Map[NN[#,bits] &,aDD];
cD=Map[NN[#,bits] &,cDD];

{mD,aD,bD,cD}                          
                           
];


InterpolationCoeff[s_,prec_]:=
Module[{aa,bb,cc,mm,nu,nuD,
interpolationData,Y,F,P,Yrule,t,
h,i,j,ai,aux,bits},

bits=Round[Log[2.,10] prec];
{mm,aa,bb,cc}=N[GaussCoefficients[s,prec],2*prec];
AppendTo[cc, 1];

interpolationData = Table[{-1+cc[[i]],Y[i]},{i,s+1}];
Yrule = {Y[i_/;i<= s]:> Y[s+1]+h Sum[(-bb[[j]]+aa[[i,j]])F[j],{j,s}]};
P[t_]=Collect[InterpolatingPolynomial[interpolationData,t],Y[_]];

aux=Table[P[cc[[i]]],{i,s}]/.Yrule//Expand;
ai=Map[NN[#,bits] &,Map[D[aux/.{h->1},#]&,Array[F,s]]]//Transpose;

nu=Table[ai[[i,j]]/bb[[j]],{i,s},{j,s}];

Rationalize[nu,0]

];



(* ::Subsection:: *)
(*End*)


End[];
EndPackage[];

