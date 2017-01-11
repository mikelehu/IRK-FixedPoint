(* ::Package:: *)

(* ::Section:: *)
(*GaussCoefficients*)


BeginPackage["GaussCoefficients`"];

NN::usage=" NN[y,m] gives m-digit precision number ";
GaussCoefficients::usage=" GaussCoefficients[s,prec] gives {m,a,b,c} s-stage Gauss method coefficients ";
InterpolationCoeff::usage=" InterpolationCoeff[s,prec] gives {nu} s-stage Gauss method's Interpolation Coefficients";

Begin["`Private`"];


(* ::Subsection:: *)
(*Functions*)


(* ::Subsubsection:: *)
(*Functions*)


NN[y_,m_]:=
(*
   y: real number.
   m: m-digit precision
*)
Module[{d,exp,naux,n},
     {d,exp}=RealDigits[y,2,m];
     naux=-exp+m;
     n=2^naux;
     Round[y*n]/n
];
SetAttributes[NN,Listable];


(* ::Text:: *)
(*Float:    32 bit -->   1 + 8 + 23  -->  eps=2^(-24)  = 6e-8.*)
(*Double: 64 bit -->  1 +11 + 52 -->   eps=2^(-53) = 1e-16.*)
(*Quadruple: 128 bit --> 1 + 15 + 112 --> eps=2^(-113) = 1e-36. *)


GaussCoefficients[s_,prec_]:=
(*

   s:     number of stages of RK method
   prec:  precision (MachinePrecision\[TildeTilde]16, 
                      QuadruplePrecision \[TildeTilde]34, ...) 

*)
Module[{highprec,i,j,
        coef,
        mantisa,
        a,b,c,m,aD,bD,cD,mD,aDD,cDD},

       highprec=200;

     coef=Flatten[NDSolve`ImplicitRungeKuttaGaussCoefficients[2*s,highprec],1];
     a=Table[coef[[i]],{i,1,s}];
     b=Table[coef[[i+s]],{i,1,s}];
     c=Table[coef[[i+2 s]],{i,1,s}];
     m=Table[a[[i,j]]/b[[j]],{i,s},{j,s}];

     mantisa=Round[Log[2.,10] prec];

     mD=Map[NN[#,mantisa] &,m];
     Do[Do[mD[[i,j]]=NN[mD[[i,j]],mantisa],{j,i-1}],{i,s}];
     Do[Do[mD[[j,i]]=1-mD[[i,j]],{j,i-1}],{i,s}];

     bD=Map[NN[#,mantisa] &,b];

     aDD= mD . DiagonalMatrix[bD];
     cDD=Table[Sum[mD[[i,j]]*bD[[j]],{j,s}],{i,s}]; 

     aD=Map[NN[#,mantisa] &,aDD];
     cD=Map[NN[#,mantisa] &,cDD];

     {mD,aD,bD,cD}                          
                           
];


InterpolationCoeff[s_,prec_]:=
(*

   s:     number of stages of RK method
   prec:  precision (MachinePrecision\[TildeTilde]16, 
                      QuadruplePrecision \[TildeTilde]34, ...) 

*)
      Module[{aa,bb,cc,mm,nu,ai,
      interpolationData,Y,F,P,Yrule,t,
      h,i,j,aux,mantisa},

      mantisa=Round[Log[2.,10] prec];
      {mm,aa,bb,cc}=N[GaussCoefficients[s,prec],2*prec];
      AppendTo[cc, 1];

      interpolationData = Table[{-1+cc[[i]],Y[i]},{i,s+1}];
      Yrule = {Y[i_/;i<= s]:> Y[s+1]+h Sum[(-bb[[j]]+aa[[i,j]])F[j],{j,s}]};
      P[t_]=Collect[InterpolatingPolynomial[interpolationData,t],Y[_]];

      aux=Table[P[cc[[i]]],{i,s}]/.Yrule//Expand;
      ai=Map[NN[#,mantisa] &,Map[D[aux/.{h->1},#]&,Array[F,s]]]//Transpose;

      nu=Table[ai[[i,j]]/bb[[j]],{i,s},{j,s}];
      Rationalize[nu,0]

];



(* ::Subsection:: *)
(*End*)


End[];
EndPackage[];

