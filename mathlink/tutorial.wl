(* ::Package:: *)

(* ::Title:: *)
(*Tutorial on FiniteFlow*)


(* ::Text:: *)
(*This is a tutorial introducing the FiniteFlow Mathematica package, which interfaces to the finiteflow C++ library.*)


(* ::Text:: *)
(*This tutorial is also beneficial in order to better understand the other publicly released examples.*)


(* ::Text:: *)
(*We obviously start by loading the package.*)


<<FiniteFlow`


(* ::Section:: *)
(*Useful wrappers*)


(* ::Text:: *)
(*The code includes some useful wrappers for common algorithms.  Among these, there are procedures for solving dense linear systems*)


FFDenseSolve[{a1 x + b1 y == c1, a2 x + b2 y == c2},{x,y}]


(* ::Text:: *)
(*sparse linear systems*)


FFSparseSolve[{a1 x + b1 y == c1, a2 x + b2 y == c2},{x,y}]


(* ::Text:: *)
(*inverting matrices*)


FFInverse[{{a,b},{c,d}}]


(* ::Text:: *)
(*and linear fits*)


(* Find c0,c1 such that
     c0/(1-x) + c1/(1+x) \[Equal] (a + b x)/(1-x^2)
   for all x, as functions of a,b. *)
FFLinearFit[{a,b},
             c0/(1-x) + c1/(1+x),
             (a + b x)/(1-x^2),
             {x},
             {c0,c1}]


(* ::Text:: *)
(*You can learn more by reading the documentation of these procedures, e.g.*)


?FFInverse


(* ::Text:: *)
(*and looking at the available Options*)


Options[FFSparseSolve]


(* ::Text:: *)
(*As an example you can speed up communication of large and sparse linear systems between Mathematica and the finiteflow library, and/or optimize the evaluation, by specifying additional options, e.g.*)


FFSparseSolve[{a1 x[1] + b1 x[2] == c1, a2 x[1] + b2 x[2] == c2},{x[1],x[2]},
              "Parameters"->{a1,b1,c1,a2,b2,c2},
              "NeededVars"->{x[1]} (* here I only care about x[1] and not x[2] *),
              "VarsPattern"->(Union[Cases[{#},_x,Infinity]]&) (*finds the list of unknowns in a expression*)
              ]


(* ::Text:: *)
(*While these wrappers can be very useful, they obviously cover only a small number of problems. We therefore strongly recommend using the approach based on dataflow graphs, which allows to define a wide variety of arbitrarily complex algorithms, and is the focus of the remainder of this tutorial.*)


(* ::Section:: *)
(*Basic usage*)


(* ::Text:: *)
(*We define a new dataflow graph and we give it an ID: mygraph.  The ID can be any Mathematica expression.  The function, if successful, returns an integer ID (which can be ignored).*)


FFNewGraph[mygraph]


(* ::Text:: *)
(*Then we define the input node, to which we assign the ID: input.  The successful definition of any node returns True.  In this example we create a graph which is a function of three variables.*)


FFGraphInputVars[mygraph,input,{x1,x2,x3}]


(* ::Text:: *)
(*Note that the names of the variables are not stored internally and they are completely irrelevant.  We might as well have used the following:*)


FFGraphInputVars[mygraph,input,{1,2,3}]


(* ::Text:: *)
(*As a shortcut, we can also define both a graph and its input node in a single command:*)


FFNewGraph[mygraph,input,{x1,x2,x3}]


(* ::Text:: *)
(*Next we want to define some algorithms in the graph, and combine them together in a more complex calculation.  Each algorithm corresponds to a node, and gets its inputs from other nodes (including the input node defined above).  Most of the commands defining algorithms have the form*)
(**)
(*  FFAlgXXX[graphname,nodename,{input1,input2,...},args...]*)
(**)
(*where XXX is the name of the algorithm.  The first argument is the ID of the graph.  The second is the ID of the new node to be created.  The third is a list of IDs corresponding to the input nodes.  Additional arguments, which are specific to the algorithm XXX,  are added later.*)


(* ::Text:: *)
(*Let us first create a node fun1 which evaluates a list of rational functions*)


FFAlgRatFunEval[mygraph,fun1,{input},
                {x1,x2,x3},
                {(1+x1)/(1+x2), (1+x2)/(1+x3), (1+x1)/(1+x3), (1+x1^2)/(1+x2^2)}]


(* ::Text:: *)
(*We won't cover all the details of all the algorithms we will use in this tutorial, but remember you can always look at the documentation of each procedure, to better understand what is going on.*)


?FFAlgRatFunEval


(* ::Text:: *)
(*Let us create another node fun2 which evaluates another list of functions*)


FFAlgRatFunEval[mygraph,fun2,{input},
                {x1,x2,x3},
                {(1+x2^2)/(1+x3^2), (1+x1^2)/(1+x3^2), (1-x1^2)/(1-x2^2), (1-x1^2)/(1-x3^2)}]


(* ::Text:: *)
(*The nodes fun1 and fun2 return a list of 4 elements each.  We now want to interpret this list as the elements of 2*2 matrices and do a matrix multiplication between them.*)


FFAlgMatMul[mygraph,matmul,{fun1,fun2},2,2,2]


(* ::Text:: *)
(*We now select this last node as the output node of our dataflow graph.*)


FFGraphOutput[mygraph,matmul]


(* ::Text:: *)
(*We can get a graphical representation of our graph with the following command: *)


FFGraphDraw[mygraph]


(* ::Text:: *)
(*Now the graph can be evaluated over finite fields.  There are a few (but not many) cases where you might want to evaluate the graph yourself, from Mathematica.  It can be done as follows:*)


FFGraphEvaluate[mygraph,{123,345,567},"PrimeNo"->0]


(* ::Text:: *)
(*The output is a list of length 4.  In  this case it represent the elements of the 2*2 matrix computed via the matrix multiplication algorithm.*)


ArrayReshape[%,{2,2}]//MatrixForm


(* ::Text:: *)
(*The option "PrimeNo" above selects the prime to be used for the evaluation from a list of hard-coded primes.  By default it is "PrimeNo"->0.  You can check its value using*)


FFPrimeNo[0]


(* ::Text:: *)
(*In most cases you will not be interested in evaluating the graph yourself, but only in reconstructing the analytic expression of its output.  This is easily done with the following command:*)


FFReconstructFunction[mygraph,{x1,x2,x3}]


(* ::Text:: *)
(*The previous command automatically performs all the needed evaluations and reconstructs the analytic formula of the output.  The evaluations and the reconstruction are automatically parallelized using multi-threading.  The "NThreads" option can be used to manually select the number of threads to be used, otherwise a choice will be made automatically depending on the hardware configuration.  We recommend specifying this option manually when using nodes in clusters or shared machines.*)


(* ::Text:: *)
(*If the graph is no longer needed, we may delete it.*)


FFDeleteGraph[mygraph]


(* ::Section:: *)
(*Linear systems and learning algorithms*)


(* ::Text:: *)
(*In this section we cover algorithms which have a learning phase.  In particular, we focus on linear systems.*)


(* ::Text:: *)
(*This example is similar to the previous one, but this time we want to compute the inverse of a matrix, and then perform the matrix multiplication.  We do not want to reconstruct the inverse analytically, but only the result of the final matrix multiplication.*)


(* ::Text:: *)
(*We create a new graph with three input variables*)


FFNewGraph[mygraph,input,{x1,x2,x3}]


(* ::Text:: *)
(*We define two analytic matrices*)


mat1 = {{(1+x2^2)/(1+x3^2), (1+x1^2)/(1+x3^2)}, {(1-x1^2)/(1-x2^2), (1-x1^2)/(1-x3^2)}};
mat2 = {{(1+x1)/(1+x2), (1+x2)/(1+x3)}, {(1+x1)/(1+x3), (1+x1^2)/(1+x2^2)}};


(* ::Text:: *)
(*We want to compute mat1.Inverse[mat2].*)


(* ::Text:: *)
(*First, we define a node computing the entries of mat1.*)


FFAlgRatFunEval[mygraph,m1,{input}, {x1,x2,x3}, Join@@mat1]


(* ::Text:: *)
(*As for the inverse of mat2, we can define it, using the Gauss-Jordan method, as the solution of the following linear system:*)
(**)
(*  mat2[[i]].{x[1],x[2]} - t[i] == 0  for i=1,2*)
(*  *)
(*with respect to the unknowns (sorted by weight) {x[1],x[2],t[1],t[2]}.  Notice that the ordering is important.*)


FFAlgDenseSolver[mygraph,invmat2,{input},{x1,x2,x3},
                 Table[mat2[[i]].{x[1],x[2]} - t[i] == 0 , {i,1,2}],
                 {x[1],x[2],t[1],t[2]}]


(* ::Text:: *)
(*We don't want the constant term in the solution (which is zero), so we exclude it with*)


FFSolverOnlyHomogeneous[mygraph,invmat2]


(* ::Text:: *)
(*A  linear solver needs to complete a learning phase, where it learns the dependent variables, independent variables, etc...   Before this is completed, the node cannot be evaluated.*)


FFGraphOutput[mygraph,invmat2];
solverlearn = FFDenseSolverLearn[mygraph,{x[1],x[2],t[1],t[2]}]


(* ::Text:: *)
(*From the lists "DepVars" and "IndepVars" we deduce that  the matrix mat is invertible.  Now we can evaluate its output.*)


(* ::Text:: *)
(*We can (optionally) check that invmat2 is indeed computing the inverse of mat2, by creating a node for mat2 and one which multiplies it with its inverse,*)


FFAlgRatFunEval[mygraph, m2, {input}, {x1,x2,x3}, Join@@mat2]

FFAlgMatMul[mygraph,check,{invmat2,m2},2,2,2]


(* ::Text:: *)
(*and then we check it by evaluating the matrix multiplication numerically*)


FFGraphOutput[mygraph,check];
ArrayReshape[FFGraphEvaluate[mygraph,{123,456,789}],{2,2}]//MatrixForm


(* ::Text:: *)
(*We now perform the matrix multiplication we said we wanted, i.e. mat1.invmat2*)


FFAlgMatMul[mygraph,matmul,{m1,invmat2},2,2,2]


(* ::Text:: *)
(*and reconstruct its output*)


FFGraphOutput[mygraph,matmul];
FFReconstructFunction[mygraph,{x1,x2,x3}]


FFDeleteGraph[mygraph]


(* ::Section:: *)
(*Subgraphs	*)


(* ::Text:: *)
(*In this section we introduce subgraph nodes.  In particular, we focus on the Laurent series expansion algorithm.*)


(* ::Text:: *)
(*We define a first graph g1 depending on three variables.  To make things simple, we just define a list of rational functions, but it could be any other arbitrarily complicated graph, whose output may not be known analytically.*)


FFNewGraph[g1,input,{x1,x2,x3}]

FFAlgRatFunEval[g1,functions,{input},{x1,x2,x3},{(1+x1 x2^2)/(x1+x1 x3^2), (1+x1^2)/(x1+x1^2 x3^2), x1 (1-x1^2)/(1-x1 x2^2), (1-x1^2)/(x1-x3^2)}]

FFGraphOutput[g1,functions]


(* ::Text:: *)
(*We now want to define a graph g2 which computes the Laurent expansion of g1 with respect to x1 (which must be the first input variable of g1) up to (and including) x^order terms, where order can be any positive, zero, or negative integer, say*)


order=3


(* ::Text:: *)
(*The coefficients of the Laurent expansion will be rational functions of x2 and x3.*)


FFNewGraph[g2,input,{x2,x3}]


(* ::Text:: *)
(*We define the node in g2 which computes the Laurent expansion, taking g1 as subgraph.*)


FFAlgLaurent[g2,laurent,{input},g1,order]


FFGraphOutput[g2,laurent]


(* ::Text:: *)
(*Laurent expansion algorithms need a learning phase, where the starting order of the expansion and other information is obtained.*)


laurentlearn = FFLaurentLearn[g2]


(* ::Text:: *)
(*Notice that we can get the general form of the expansion, in terms of generic coefficients,*)


coefficients = coeff/@Range[FFNParsOut[g2,laurent]]


(* ::Text:: *)
(*based on the information we already learned, and without doing any reconstruction.*)


FFLaurentSol[coefficients,x1,laurentlearn]


(* ::Text:: *)
(*The coefficients of the Laurent expansion can be computed numerically and reconstructed.  *)


reconstructed=FFReconstructFunction[g2,{x2,x3}]


(* ::Text:: *)
(*From the reconstructed coefficients we can explicitly build the analytic Laurent expansion using the FFLaurentSol utility*)


FFLaurentSol[reconstructed,x1,laurentlearn]


(* ::Text:: *)
(*Notice that, as mentioned, the coefficients of the Laurent expansion computed by g2 can be obtained without knowing the analytic expression of the output of g1.*)
