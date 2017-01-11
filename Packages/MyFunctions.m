(* ::Package:: *)

(* ::Section:: *)
(*MyFunctions*)


BeginPackage["MyFunctons`"];

(*******MyFunctions-1**********)

ErrorPosition::usage=" ErrorPosition[y1,y2,prec] gives Norm[q1-q2] 
                       where q1,q2 are position vectors of corresponding solution";

Chdata::usage=" Chdata[u,parameters] gives unew coordinates values to 
                get zero linear momentum";

Pert::usage=" Pert[e0,k] gives random perturbations using k integer.";
FunT12::usage= " ";

(*******MyFunctions-2**********)

FunEnergy::usage=" Gives Mean-Energy-Error and Desviation-Energy-Error of input integrations"; 
FunAllEnergy::usage=" "; 
FunError::usage=" Gives the position error respect a reference solution"; 
FunErrorNEW::usage=" Gives the position error respect a reference solution"; 
FunHamErrDif::usage=" Gives "; 
FunEstimation::usage=" "; 

FunMomentum::usage=" Gives Mean-Angular-Momentum-Error and Desviation-Angular-Momentum-Error of input integrations"; 

(******* Fortran esperiments *******)

FunEnergyFortran::usage=" Gives Mean-Energy-Error and Desviation-Energy-Error of input integrations"


Begin["`Private`"];




(* ::Subsection::Closed:: *)
(*MyFunctions-1*)


ErrorPosition[y1_,y2_ ,prec_:100]:=
 (*
   y1 = (q_1, p_1, eq_1, ep_1);
   y2 = (q_2, p_2, eq_2, ep_2);
 *)

Module[{q1,eq1,q2,eq2,neq,d},

        neq=(Length[y1]-1)/2;
        d=neq/2;

        q1=SetPrecision[y1[[2;;d+1]],prec];
        eq1=SetPrecision[y1[[2d+2;;3d+1]],prec];
        q2=SetPrecision[y2[[2;;d+1]],prec];
        eq2=SetPrecision[y2[[2d+2;;3d+1]],prec];

        Norm[(q1+eq1)-(q2+eq2)]
      ];




Chdata[u_,parameters_List]:=
  (* 
     N-Body problem.
     u = (q,p);
     parameters= Gm_i, i=1,N
   *)
Module[{dim,neq,d,nbody,Gm,MM,mq,mv,i,i1,i2,id,qi,vi,qnew,vnew},

        dim=3;
        neq=Length[u];
        d=neq/2;
        nbody=Length[u]/(2*dim);
        Gm=parameters;
        MM=Sum[Gm[[i]],{i,nbody}];
        mq=Array[0&,dim];
        mv=Array[0&,dim];

        Do
        [  i1=(i-1)*dim;
           i2=d+i1;
           qi=Table[u[[i1+id]],{id,1,dim}];
           vi=Table[u[[i2+id]],{id,1,dim}];
           mq=mq+Gm[[i]]*qi;
           mv=mv+Gm[[i]]*vi,
        {i,nbody}
        ];

        mq=mq/MM;
        mv=mv/MM;

        qnew={};
        vnew={};

        Do
        [  i1=(i-1)*dim;
           i2=d+i1;
           AppendTo[qnew,Table[u[[i1+id]]-mq[[id]],{id,1,dim}]];
           AppendTo[vnew,Table[u[[i2+id]]-mv[[id]],{id,1,dim}]],
        {i,nbody}
        ];

        Flatten[Join[qnew,vnew]]
];



Pert[e0_,k_]:={e0*(k*RandomReal[{-1,1}])};


FunT12[alpha_,t_]:=alpha*Sqrt[t];



(* ::Subsection:: *)
(*MyFunctions-2*)


FunEnergy[outAll_, nstat_,nout_,HAM_,parameters_
             ,neq_,prec_,step_Integer]:=
(*
      outALL: file with solution of all integrations.
      nstat: number of integrations.
      nout: number of mesh points of each integration.
      HAM: name of Hamiltonian function.
      parameters: parameters for Hamiltonian.
      neq: number of equations.
      prec: precision.
      step:  
*)


Module[{i,outA,Ham0,Hamerr,
        SHamerr=Array[0 &,Floor[nout/step]],
        DHamerr=Array[0 &,Floor[nout/step]],
        MeanE,DesvE},

For[i=1, i<= nstat,i++,
  outA=outAll[[i]];
  Ham0 = HAM[neq,First[outA],parameters,prec];
  Hamerr = Map[HAM[neq,#,parameters,prec]/Ham0-1&,outA];
  Hamerr=  Map[First,Partition[Hamerr,step]];
  SHamerr=SHamerr+Hamerr;
  DHamerr=DHamerr+Hamerr^2;
 
];
MeanE=SHamerr/nstat;
DesvE=Sqrt[DHamerr/nstat-MeanE^2];
{MeanE,DesvE}
]



FunAllEnergy[outAll_, nstat_,nout_,HAM_,parameters_
             ,neq_,prec_,step_Integer]:=
Module[{i, outA,
        Ham0, Hamerr,AllHamerr={}
},
For[i=1, i<= nstat,i++,
  outA=outAll[[i]];
  Ham0 = HAM[neq,First[outA],parameters,prec];
  Hamerr = Map[HAM[neq,#,parameters,prec]/Ham0-1&,outA];
  Hamerr=  Map[First,Partition[Hamerr,step]];
  AppendTo[AllHamerr,Hamerr];
];
AllHamerr
]


FunError[outRef_, outAll_,ERR_, nstat_,nout_,
         neq_,prec_,step_Integer]:=
(*
   outRef: name of the file with reference solutions.
   outAll: name of the file with solutions to evaluate the error.
   ERR: name of the mathematica function to evaluate the error.
   nstat: number of integrations.
   nout: number of mesh points of each integration.
   neq:
   prec:
   step:
*)

Module[{i,Ref,outA,
        Error,
        SError=Array[0 &,Floor[nout/step]]},
For[i=1, i<= nstat,i++,
  Ref=outRef[[i]];
  outA=outAll[[i]];
  Error =Chop[ MapThread[ERR[#1,#2,prec]&,{Ref, outA}],10^-100];
  Error= Map[First,Partition[Error,step]];
  SError=SError+Abs[Error];
];
SError/nstat
]



FunErrorNEW[outRef_, outAll_,ERR_, nstat_,nout_,
         neq_,prec_,step_Integer]:=
(*
   outRef: name of the file with reference solutions.
   outAll: name of the file with solutions to evaluate the error.
   ERR: name of the mathematica function to evaluate the error.
   nstat: number of integrations.
   nout: number of mesh points of each integration.
   neq:
   prec:
   step:
*)

Module[{i,Ref,outA,
        Error,
        SError=Array[0 &,Floor[nout/step]],
        DError=Array[0 &,Floor[nout/step]],
        MeanE,DesvE},
        
For[i=1, i<= nstat,i++,
  Ref=outRef[[i]];
  outA=outAll[[i]];
  Error =Chop[ MapThread[ERR[#1,#2,prec]&,{Ref, outA}],10^-100];
  Error= Map[First,Partition[Error,step]];
  SError=SError+Abs[Error];
  DError=DError+Error^2;
  
];
MeanE=SError/nstat;
DesvE=Sqrt[DError/nstat-MeanE^2];
{MeanE,DesvE}
]


FunHamErrDif[outAll_, nstat_,nout_,HAM_,parameters_,
             neq_,prec_,step_Integer]:=
Module[{i,Ham0, Hamerr,Hamerrdif,Error,outA},
Hamerrdif={};
For[i=1, i<= nstat,i++,
  outA=outAll[[i]];
  Ham0 = HAM[neq,First[outA],parameters,prec];
  Hamerr=Map[HAM[neq,#,parameters,prec]/Ham0-1&,outA];
  Hamerr=Map[First,Partition[Hamerr,step]];
  Hamerrdif = Join[Hamerrdif,Drop[Hamerr,1]-Drop[Hamerr,-1]];
];
Hamerrdif
]


FunEstimation[outRef_,outAll_,ERR_,nstat_,nout_,
              neq_,prec_,step_Integer]:=
(*
   outRef: name of the file with reference solutions.
   outAll: name of the file with estimations of the integrations.
   ERR: name of the mathematica function to evaluate the error.
   nstat:
   nout:
   neq:
   prec:
   step:
*)

Module[{i,outA,Ref,
        est,estQ,meanEst,desvEst,Error,
        Qty,meanQty,desvQty,
        Sout=Array[0 &,Floor[nout/step]],
        S2out=Array[0 &,Floor[nout/step]],
        SEstQty=Array[0 &,Floor[nout/step]],
        S2EstQty=Array[0 &,Floor[nout/step]]},

For[i=1, i<= nstat,i++,
  outA=outAll[[i]];
  Ref=outRef[[i]];
  est=Take[outA, All,{2neq+2,3neq+1}];
  estQ=Map[Norm,Chop[Take[est,All,{1,neq/2}],10^-100]];
  estQ=Map[First,Partition[estQ,step]];
  Sout=Sout+estQ;
  S2out=S2out+estQ^2;
  Error =Chop[ MapThread[ERR[#1,#2,prec]&,{Ref, outA}],10^-100];
  Error= Map[First,Partition[Error,step]];
  Qty=Log[10,Quiet[estQ/Error]/.ComplexInfinity->1/.Indeterminate->1//.0->1];
  SEstQty=SEstQty+Qty;
  S2EstQty=S2EstQty+Qty^2;
];

meanEst=Sout/nstat;
desvEst=Sqrt[S2out/nstat-meanEst^2];
meanQty=SEstQty/nstat;
desvQty=Sqrt[S2EstQty/nstat-meanQty^2];
{meanEst,desvEst,meanQty,desvQty}
];


FunMomentum[outAll_, nstat_,nout_,MOM_,parameters_
             ,neq_,prec_,step_Integer]:=
(*
      outALL: file with solution of all integrations.
      nstat: number of integrations.
      nout: number of mesh points of each integration.
      MOM: name of Angular-Momentum function.
      parameters: parameters for Hamiltonian.
      neq: number of equations.
      prec: precision.
      step:  
*)

Module[{i,outA,L0,Lerr,
        SLerr=Array[0 &,Floor[nout/step]],
        DLerr=Array[0 &,Floor[nout/step]],
        MeanL,DesvL},

For[i=1, i<= nstat,i++,
  outA=outAll[[i]];
  L0 = MOM[neq,First[outA],parameters,prec];
  Lerr = Map[MOM[neq,#,parameters,prec]/L0-1&,outA];
  Lerr=  Map[First,Partition[Lerr,step]];
  SLerr=SLerr+Lerr;
  DLerr=DLerr+Lerr^2;
 
];
MeanL=SLerr/nstat;
DesvL=Sqrt[DLerr/nstat-MeanL^2];
{MeanL,DesvL}
]


(* ::Subsection::Closed:: *)
(*MyFunctions: Fortran*)


FunEnergyFortran[outAll_, nstat_,nout_,HAM_,parameters_
                 ,neq_,prec_,step_Integer, LimitErr_]:=
(*
      outALL: file with solution of all integrations.
      nstat: number of integrations.
      nout: number of mesh points of each integration.
      HAM: name of Hamiltonian function.
      parameters: parameters for Hamiltonian.
      neq: number of equations.
      prec: precision.
      step:  
*)

Module[{i,k=0,outA,Ham0,Hamerr,
        SHamerr=Array[0 &,Floor[nout/step]],
        DHamerr=Array[0 &,Floor[nout/step]],
        MeanE,DesvE},

For[i=1, i<= nstat,i++,
  outA=outAll[[i]];
  Ham0 = HAM[neq,First[outA],parameters,prec];
  Hamerr = Map[HAM[neq,#,parameters,prec]/Ham0-1&,outA];
  Hamerr=  Map[First,Partition[Hamerr,step]];
  If[Max[Abs[Hamerr]]<LimitErr,
       SHamerr=SHamerr+Hamerr;
       DHamerr=DHamerr+Hamerr^2,
     k++];
 
];
Print[k];
MeanE=SHamerr/(nstat-k);
DesvE=Sqrt[DHamerr/(nstat-k)-MeanE^2];
{MeanE,DesvE}
]


FunAllEnergy[outAll_, nstat_,nout_,HAM_,parameters_
             ,neq_,prec_,step_Integer]:=
Module[{i, outA,
        Ham0, Hamerr,AllHamerr={}
},
For[i=1, i<= nstat,i++,
  outA=outAll[[i]];
  Ham0 = HAM[neq,First[outA],parameters,prec];
  Hamerr = Map[HAM[neq,#,parameters,prec]/Ham0-1&,outA];
  Hamerr=  Map[First,Partition[Hamerr,step]];
  AppendTo[AllHamerr,Hamerr];
];
AllHamerr
]


FunError[outRef_, outAll_,ERR_, nstat_,nout_,
         neq_,prec_,step_Integer]:=
(*
   outRef: name of the file with reference solutions.
   outAll: name of the file with solutions to evaluate the error.
   ERR: name of the mathematica function to evaluate the error.
   nstat: number of integrations.
   nout: number of mesh points of each integration.
   neq:
   prec:
   step:
*)

Module[{i,Ref,outA,
        Error,
        SError=Array[0 &,Floor[nout/step]]},
For[i=1, i<= nstat,i++,
  Ref=outRef[[i]];
  outA=outAll[[i]];
  Error =Chop[ MapThread[ERR[#1,#2,prec]&,{Ref, outA}],10^-100];
  Error= Map[First,Partition[Error,step]];
  SError=SError+Abs[Error];
];
SError/nstat
]


FunHamErrDif[outAll_, nstat_,nout_,HAM_,parameters_,
             neq_,prec_,step_Integer]:=
Module[{i,Ham0, Hamerr,Hamerrdif,Error,outA},
Hamerrdif={};
For[i=1, i<= nstat,i++,
  outA=outAll[[i]];
  Ham0 = HAM[neq,First[outA],parameters,prec];
  Hamerr=Map[HAM[neq,#,parameters,prec]/Ham0-1&,outA];
  Hamerr=Map[First,Partition[Hamerr,step]];
  Hamerrdif = Join[Hamerrdif,Drop[Hamerr,1]-Drop[Hamerr,-1]];
];
Hamerrdif
]


FunEstimation[outRef_,outAll_,ERR_,nstat_,nout_,
              neq_,prec_,step_Integer]:=
(*
   outRef: name of the file with reference solutions.
   outAll: name of the file with estimations of the integrations.
   ERR: name of the mathematica function to evaluate the error.
   nstat:
   nout:
   neq:
   prec:
   step:
*)

Module[{i,outA,Ref,
        est,estQ,meanEst,desvEst,Error,
        Qty,meanQty,desvQty,
        Sout=Array[0 &,Floor[nout/step]],
        S2out=Array[0 &,Floor[nout/step]],
        SEstQty=Array[0 &,Floor[nout/step]],
        S2EstQty=Array[0 &,Floor[nout/step]]},

For[i=1, i<= nstat,i++,
  outA=outAll[[i]];
  Ref=outRef[[i]];
  est=Take[outA, All,{2neq+2,3neq+1}];
  estQ=Map[Norm,Chop[Take[est,All,{1,neq/2}],10^-100]];
  estQ=Map[First,Partition[estQ,step]];
  Sout=Sout+estQ;
  S2out=S2out+estQ^2;
  Error =Chop[ MapThread[ERR[#1,#2,prec]&,{Ref, outA}],10^-100];
  Error= Map[First,Partition[Error,step]];
  Qty=Log[10,Quiet[estQ/Error]/.ComplexInfinity->1/.Indeterminate->1//.0->1];
  SEstQty=SEstQty+Qty;
  S2EstQty=S2EstQty+Qty^2;
];

meanEst=Sout/nstat;
desvEst=Sqrt[S2out/nstat-meanEst^2];
meanQty=SEstQty/nstat;
desvQty=Sqrt[S2EstQty/nstat-meanQty^2];
{meanEst,desvEst,meanQty,desvQty}
];


FunMomentum[outAll_, nstat_,nout_,MOM_,parameters_
             ,neq_,prec_,step_Integer]:=
(*
      outALL: file with solution of all integrations.
      nstat: number of integrations.
      nout: number of mesh points of each integration.
      MOM: name of Angular-Momentum function.
      parameters: parameters for Hamiltonian.
      neq: number of equations.
      prec: precision.
      step:  
*)

Module[{i,outA,L0,Lerr,
        SLerr=Array[0 &,Floor[nout/step]],
        DLerr=Array[0 &,Floor[nout/step]],
        MeanL,DesvL},

For[i=1, i<= nstat,i++,
  outA=outAll[[i]];
  L0 = MOM[neq,First[outA],parameters,prec];
  Lerr = Map[MOM[neq,#,parameters,prec]/L0-1&,outA];
  Lerr=  Map[First,Partition[Lerr,step]];
  SLerr=SLerr+Lerr;
  DLerr=DLerr+Lerr^2;
 
];
MeanL=SLerr/nstat;
DesvL=Sqrt[DLerr/nstat-MeanL^2];
{MeanL,DesvL}
]


(* ::Subsection:: *)
(*End*)


End[ ];
EndPackage[];



