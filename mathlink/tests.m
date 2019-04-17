(* ::Package:: *)

<<FiniteFlow`


Print["Test multivariate reconstruction"];
FFNewGraph[g];
vars = {x1,x2,x3};
funcs = Join@@Table[{(1 + x1^60 + x2^60 + x3^60)(1+x1+x2),
               ii + 2^(8 ii) x1 + 2^64 x2,
               2^(8 ii) + x1 + x2 + x3},{ii,10}];
FFGraphInputVars[g,in,vars];
FFAlgRatFunEval[g,rf,{in},vars,funcs];
FFGraphOutput[g,rf];
res = FFReconstructFunction[g,{x1,x2,x3}];
If[!AllTrue[res-funcs /. {x1->12/19,x2->30/17,x3->23/41},#==0&],
    Print["- Test FAILED!"],
    Print["- Test passed"];
  ];
FFDeleteGraph[g];






