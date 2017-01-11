(* ::Package:: *)

(* ::Section:: *)
(*N-Body Problem *)


BeginPackage["NBodyProblem`"];


NBodyODE::usage=" ";
NBodyODEc::usage=" ";
NBodyODEd::usage=" ";
NBodyODEfunc::usage=" ";

NBodyHam::usage="";
NBodyMom::usage="";

NBodyHam3::usage=" Fortran";

Begin["`Private`"];


(* ::Subsection:: *)
(*Functions*)


(* ::Subsubsection::Closed:: *)
(*Odefun*)


NBodyODE[t_,u_,parameters_List]:=
Module[{evalODE,dim,neq,d,nbody,y,r3,Gm,f,f2,i,i1,i2,j,j1,j2,id,
             qi,qj,qij,vi},
dim=3;
neq=Length[u];
d=neq/2;
nbody=Length[u]/(2*dim);
f=Array[0 &,neq];
f2=Array[0 &,d];
Gm =Drop[parameters,-1];

If[parameters[[-1]]==1, 
(*fq*)
Do[ i1=(i-1)*dim;
      i2=d+i1;
      vi=Table[u[[i2+id]],{id,1,dim}];
      Do[f2[[i1+id]]=vi[[id]],{id,1,dim}],     
{i,nbody}];
f2,

If[parameters[[-1]]==2,
(*fp*)
Do[ i1=(i-1)*dim;
      qi=Table[u[[i1+id]],{id,1,dim}];
      Do [    j1=(j-1)*dim;
                 qj=Table[u[[j1+id]],{id,1,dim}];
                 qij=qi-qj;
                 r3=Norm[qij]^3;
                 Do[ f2[[i1+id]]=f2[[i1+id]]-Gm[[j]]*qij[[id]]/r3,{id,1,dim}];
                 Do[ f2[[j1+id]]=f2[[j1+id]]+Gm[[i]]*qij[[id]]/r3,{id,1,dim}],
 {j,i+1,nbody}],
{i,nbody}];
f2,

(*else fq,fp*)
Do[ i1=(i-1)*dim;
      i2=d+i1;
      qi=Table[u[[i1+id]],{id,1,dim}];
      vi=Table[u[[i2+id]],{id,1,dim}];
      Table[f[[i1+id]]=vi[[id]],{id,1,dim}];
      Do [    j1=(j-1)*dim;
                j2=d+j1;
                 qj=Table[u[[j1+id]],{id,1,dim}];
                 qij=qi-qj;
                 r3=Norm[qij]^3;
                 Do[ f[[i2+id]]=f[[i2+id]]-Gm[[j]]*qij[[id]]/r3,{id,1,dim}];
                 Do[ f[[j2+id]]=f[[j2+id]]+Gm[[i]]*qij[[id]]/r3,{id,1,dim}],
 {j,i+1,nbody}],
{i,nbody}];
f
]]

];

NBodyODEc[t_,u_,parameters_List]:=Apply[NBodyODEfunc,Flatten[{t,u,parameters}]];

Quadprec=Abs[Log[10.,2^-113]];
NBodyODEd[t_,u_,parameters_List]:=SetPrecision[Apply[NBodyODEfunc,Flatten[N[{t,u,parameters}]]],Quadprec];

NBodyODEfunc=Compile[{t,
q11,q12,q13,q21,q22,q23,q31,q32,q33,q41,q42,q43,q51,q52,q53,q61,q62,q63,q71,q72,q73,q81,q82,q83,q91,q92,q93,q101,q102,q103,
p11,p12,p13,p21,p22,p23,p31,p32,p33,p41,p42,p43,p51,p52,p53,p61,p62,p63,p71,p72,p73,p81,p82,p83,p91,p92,p93,p101,p102,p103,
Gma1,Gma2,Gma3,Gma4,Gma5,Gma6,Gma7,Gma8,Gma9,Gma10,par1},

Module[{dim,neq,d,nbody,y,r3,Gm,f,f2,i,i1,i2,j,j1,j2,id,
             qi,qj,qij,vi,
             u},

u={q11,q12,q13,q21,q22,q23,q31,q32,q33,q41,q42,q43,q51,q52,q53,q61,q62,q63,q71,q72,q73,q81,q82,q83,q91,q92,q93,q101,q102,q103,
p11,p12,p13,p21,p22,p23,p31,p32,p33,p41,p42,p43,p51,p52,p53,p61,p62,p63,p71,p72,p73,p81,p82,p83,p91,p92,p93,p101,p102,p103};

dim=3;
neq=60;
d=30;(*neq/2       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*)
nbody=neq/(2*dim);
f=Array[N[0] &,neq];
f2=Array[N[0] &,d];
Gm ={Gma1,Gma2,Gma3,Gma4,Gma5,Gma6,Gma7,Gma8,Gma9,Gma10};

If[par1==1, 
(*fq*)
Do[ i1=(i-1)*dim;
      i2=d+i1;
      vi=Table[u[[i2+id]],{id,1,dim}];
      Do[f2[[i1+id]]=vi[[id]],{id,1,dim}],     
{i,nbody}];
f2,

If[par1==2,
(*fp*)
Do[ i1=(i-1)*dim;
      qi=Table[u[[i1+id]],{id,1,dim}];
      Do [    j1=(j-1)*dim;
                 qj=Table[u[[j1+id]],{id,1,dim}];
                 qij=qi-qj;
                 r3=Norm[qij]^3;
                 Do[ f2[[i1+id]]=f2[[i1+id]]-Gm[[j]]*qij[[id]]/r3,{id,1,dim}];
                 Do[ f2[[j1+id]]=f2[[j1+id]]+Gm[[i]]*qij[[id]]/r3,{id,1,dim}],
 {j,i+1,nbody}],
{i,nbody}];
f2,

(*else fq,fp*)
Do[ i1=(i-1)*dim;
      i2=d+i1;
      qi=Table[u[[i1+id]],{id,1,dim}];
      vi=Table[u[[i2+id]],{id,1,dim}];
      Do[f[[i1+id]]=vi[[id]],{id,1,dim}];
      Do [    j1=(j-1)*dim;
                j2=d+j1;
                 qj=Table[u[[j1+id]],{id,1,dim}];
                 qij=qi-qj;
                 r3=Norm[qij]^3;
                 Do[ f[[i2+id]]=f[[i2+id]]-Gm[[j]]*qij[[id]]/r3,{id,1,dim}];
                 Do[ f[[j2+id]]=f[[j2+id]]+Gm[[i]]*qij[[id]]/r3,{id,1,dim}],
 {j,i+1,nbody}],
{i,nbody}];
f
]]
]
];



(* ::Subsubsection::Closed:: *)
(*Hamiltonian C*)


NBodyHam[neq_,out_,parameters_List, doi_:100]:=
Module[{dim,d,nbody,Gm,uu,H,P,r,i,i1,i2,j,j1,j2,id,
              qi,qj,qij,vi,
              d1,d2,d3,d4},
dim=3;
d=neq/2;
nbody=neq/(2*dim);

d1=2;
d2=nbody*(2*dim)+1;
d3=d2+1;
d4=d3+nbody*(2*dim)-1;

uu=SetPrecision[out[[d1;;d2]],doi]+SetPrecision[out[[d3;;d4]],doi];

Gm =SetPrecision[parameters,doi];
H=0;
P=0;
Do[ i1=(i-1)*dim;
      i2=d+i1;
      qi=Table[uu[[i1+id]],{id,1,dim}];
      vi=Table[uu[[i2+id]],{id,1,dim}];
      H=H+Gm[[i]]*(vi.vi);
      Do [    j1=(j-1)*dim;
                j2=d+j1;
                qj=Table[uu[[j1+id]],{id,1,dim}];
                qij=qi-qj;
                r=Norm[qij];
                P=P+Gm[[i]]*Gm[[j]]/r,
 {j,i+1,nbody}],
{i,nbody}];
H/2-P
];


(* ::Subsubsection:: *)
(*Hamiltonian Fortran*)


NBodyHam3[neq_,out_,parameters_List, doi_:100]:=
Module[{dim,d,nbody,Gm,uu,H,P,r,i,i1,i2,j,j1,j2,id,
              qi,qj,qij,vi,
              d1,d2,d3,d4},
dim=3;
d=neq/2;
nbody=neq/(2*dim);

d1=2;
d2=nbody*(2*dim)+1;
d3=d2+1;
d4=d3+nbody*(2*dim)-1;

(*uu=SetPrecision[out[[d1;;d2]],doi]+SetPrecision[out[[d3;;d4]],doi];*)
uu=SetPrecision[out[[d1;;d2]],doi];

Gm =SetPrecision[parameters,doi];
H=0;
P=0;
Do[ i1=(i-1)*dim;
      i2=d+i1;
      qi=Table[uu[[i1+id]],{id,1,dim}];
      vi=Table[uu[[i2+id]],{id,1,dim}];
      H=H+Gm[[i]]*(vi.vi);
      Do [    j1=(j-1)*dim;
                j2=d+j1;
                qj=Table[uu[[j1+id]],{id,1,dim}];
                qij=qi-qj;
                r=Norm[qij];
                P=P+Gm[[i]]*Gm[[j]]/r,
 {j,i+1,nbody}],
{i,nbody}];
H/2-P
];


(* ::Subsection::Closed:: *)
(*Angular Momentum*)


NBodyMom[neq_,out_,parameters_List, doi_:100]:=
Module[{dim,d,nbody,Gm,qv,
        d1,d2,d3,d4},
dim=3;
d=neq/2;
nbody=neq/(2*dim);

d1=2;
d2=nbody*(2*dim)+1;
d3=d2+1;
d4=d3+nbody*(2*dim)-1;

Gm =SetPrecision[parameters,doi];
qv=Partition[SetPrecision[out[[d1;;d2]],doi]+SetPrecision[out[[d3;;d4]],doi],dim];

Sum[Gm[[i]]*(qv[[i]]\[Cross]qv[[nbody+i]]),{i,nbody}]
];


(* ::Subsection:: *)
(*End*)


End[];
EndPackage[];



