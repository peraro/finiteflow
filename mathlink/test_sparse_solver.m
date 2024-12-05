(* ::Package:: *)

<<FiniteFlow`


(* solve a random system with some zeroes in the solution *)


sys={{y zb[0]+y zc[0]+y zd[0]+y ze[0]==0,y zb[5]+y zc[5]+y zd[5]+y ze[5]==0,y zc[0]+y zd[0]+y zf[0]==0,za[0]+y zb[0]+y zb[4]+y zc[0]+y zc[4]+y zc[5]+y zd[4]+y zd[5]+y ze[4]+y zf[5]==0,za[5]+y zb[5]+y zc[5]==0,y zc[4]+y zd[4]+y zf[4]==0,za[4]+y zb[4]+y zc[4]==0,y zb[0]+y ze[0]+y zf[0]==0,za[0]+y zb[0]+y zb[3]+y zb[5]+y zc[0]+y zc[3]+y zd[3]+y ze[3]+y ze[5]+y zf[5]==0,za[5]+y zb[5]+y zc[5]==0,za[0]+y zb[0]+y zb[4]+y zc[0]+y zc[3]+y zd[3]+y ze[4]+y zf[3]+y zf[4]==0,za[3]+za[4]+za[5]+y zb[3]+y zb[4]+y zb[5]+y zc[3]+y zc[4]+y zc[5]==0,za[4]+y zb[4]+y zc[4]==0,y zb[3]+y ze[3]+y zf[3]==0,za[3]+y zb[3]+y zc[3]==0,za[3]+y zb[3]+y zc[3]==0,y zb[0]+y ze[0]+y zf[0]==0,za[0]+y zb[2]+y zb[5]+y zc[2]+y zd[0]+y zd[2]+y ze[0]+y ze[2]+y ze[5]+y zf[5]==0,za[5]+y zd[5]+y ze[5]==0,za[0]+y zb[0]+y zb[4]+y zc[2]+y zd[0]+y zd[2]+y ze[4]+y zf[0]+y zf[2]+y zf[4]==0,za[0]+za[2]+za[4]+za[5]+y zb[2]+y zb[5]+y zc[2]+y zd[4]+y zd[5]+y ze[4]+y zf[5]==0,za[5]==0,za[4]+y zb[4]+y zd[4]+y zf[4]==0,za[4]==0,y zb[0]+y zb[2]+y zb[3]+y ze[0]+y ze[2]+y ze[3]+y zf[0]+y zf[2]+y zf[3]==0,za[0]+za[2]+za[3]+y zb[2]+y zb[5]+y zc[2]+y zd[3]+y ze[3]+y ze[5]+y zf[5]==0,za[5]==0,za[0]+za[2]+za[3]+y zb[2]+y zb[3]+y zb[4]+y zc[2]+y zd[3]+y ze[4]+y zf[3]+y zf[4]==0,za[3]+za[4]+za[5]==0,za[4]==0,y zb[3]+y ze[3]+y zf[3]==0,za[3]==0,za[3]==0,y zb[2]+y ze[2]+y zf[2]==0,za[2]+y zd[2]+y ze[2]==0,za[2]+y zb[2]+y zd[2]+y zf[2]==0,za[2]==0,y zb[2]+y ze[2]+y zf[2]==0,za[2]==0,za[2]==0,y zc[0]+y zd[0]+y zf[0]==0,za[0]+y zb[1]+y zc[1]+y zc[5]+y zd[0]+y zd[1]+y zd[5]+y ze[0]+y ze[1]+y zf[5]==0,za[5]+y zd[5]+y ze[5]==0,y zc[0]+y zc[1]+y zc[4]+y zd[0]+y zd[1]+y zd[4]+y zf[0]+y zf[1]+y zf[4]==0,za[0]+za[1]+za[4]+y zb[1]+y zc[1]+y zc[5]+y zd[4]+y zd[5]+y ze[4]+y zf[5]==0,za[5]==0,y zc[4]+y zd[4]+y zf[4]==0,za[4]==0,za[0]+y zb[1]+y zc[0]+y zc[3]+y zd[3]+y ze[0]+y ze[1]+y zf[0]+y zf[1]+y zf[3]==0,za[0]+za[1]+za[3]+za[5]+y zb[1]+y zc[1]+y zc[5]+y zd[3]+y ze[3]+y ze[5]+y zf[5]==0,za[5]==0,za[0]+za[1]+za[4]+y zb[1]+y zc[1]+y zc[3]+y zc[4]+y zd[3]+y ze[4]+y zf[3]+y zf[4]==0,za[3]+za[4]+za[5]==0,za[4]==0,za[3]+y zc[3]+y ze[3]+y zf[3]==0,za[3]==0,za[3]==0,za[0]+y zb[1]+y zc[2]+y zd[0]+y zd[2]+y ze[0]+y ze[1]+y zf[1]+y zf[2]==0,za[1]+za[2]+za[5]+y zd[1]+y zd[2]+y zd[5]+y ze[1]+y ze[2]+y ze[5]==0,za[0]+za[1]+za[4]+y zb[1]+y zc[2]+y zd[1]+y zd[2]+y zd[4]+y ze[4]+y zf[1]+y zf[2]==0,za[1]+za[2]+za[5]==0,za[4]==0,za[0]+za[2]+za[3]+y zb[1]+y zc[2]+y zd[3]+y ze[1]+y ze[2]+y ze[3]+y zf[1]+y zf[2]==0,za[1]+za[2]+za[5]==0,za[1]+za[2]+za[3]+za[4]==0,za[3]==0,za[2]+y zd[2]+y ze[2]==0,za[2]==0,za[2]==0,y zc[1]+y zd[1]+y zf[1]==0,za[1]+y zd[1]+y ze[1]==0,y zc[1]+y zd[1]+y zf[1]==0,za[1]==0,za[1]+y zc[1]+y ze[1]+y zf[1]==0,za[1]==0,za[1]==0,za[1]+y zd[1]+y ze[1]==0,za[1]==0,za[1]==0},{za[1],za[2],za[3],za[4],za[5],za[0],zb[1],zb[2],zb[3],zb[4],zb[5],zb[0],zc[1],zc[2],zc[3],zc[4],zc[5],zc[0],zd[1],zd[2],zd[3],zd[4],zd[5],zd[0],ze[1],ze[2],ze[3],ze[4],ze[5],ze[0],zf[1],zf[2],zf[3],zf[4],zf[5],zf[0]}};


sol = FFSparseSolve@@sys;


If[!TrueQ[sol[[0]]==List]||!TrueQ[Union[Expand[sys[[1]] /. Dispatch[sol]]]=={True}],
  Print["ERROR: FFSparseSolve test failed"];
  Quit[];
];


(*  *)


check=FFSparseSolve[
system={
x + 2 y + 3 z == a,
7 x + 3 y - 4 z == b,
x + 3 y + 13 z == c
},
vars={x,y,z}]


(* *)


TestMaxCol[maxcol_,backsubst_,keepfullout_]:=Module[
{tcheck},
tcheck[$Failed]:=(
  Print["TestMaxCol failed with inputs ",{maxcol,backsubst,keepfullout}];
  Quit[];
);
tcheck[a_]:=a;
FFNewGraph["g","in",params={a,b,c}] // tcheck;
FFAlgSparseSolver["g","sol",{"in"},params,
  system,
  vars,
  "NeededVars"->vars[[2;;]]
] // tcheck;
FFSolverSparseOutputWithMaxCol[
  "g","sol",maxcol,
  "BackSubstitution"->backsubst,
  "KeepFullOutput"->keepfullout
] // tcheck;
FFGraphOutput["g","sol"] // tcheck;
learn = FFSparseSolverLearn["g",{x,y,z}] // tcheck;
FFSparseSolverMarkAndSweepEqs["g","sol"] // tcheck;
FFSparseSolverDeleteUnneededEqs["g","sol"] // tcheck;
partsol = FFSparseSolverSol[
  FFReconstructFunction["g",params],
  learn
] // tcheck;
If[!TrueQ[Union[partsol[[;;,1]]-partsol[[;;,2]]/.check//Expand]=={0}],
  check[$Failed];
];
FFDeleteGraph[graph]//check;
]


TestMaxCol[1,False,True];
TestMaxCol[2,False,True];
TestMaxCol[1,True,True];
TestMaxCol[2,True,True];
TestMaxCol[1,True,False];
TestMaxCol[2,True,False];
TestMaxCol[1,False,False];
TestMaxCol[2,False,False];


Print["All tests passed!"];


Quit[];
