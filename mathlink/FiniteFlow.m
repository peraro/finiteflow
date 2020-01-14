(* ::Package:: *)

BeginPackage["FiniteFlow`"]


FFInt64Max::usage = "Maximum 64-bit integer."

FFLoadLib::usage = "FFLoadLib[] (re)loads the finiteflow library."

FFBadShift::usage = "Returned by the reconstruction routines when there is a bad shift."
FFImpossible::usage = "Returned by the solve routines when the system is impossible."
FFSingularMatrix::usage = "Returned by FFInverse when the input matrix is singular."

FFCoefficientRules::usage = "FFCoefficientRules[expr,vars] is an alternative to the builtin CoefficientRules[expr,vars] which uses CoefficientList internally."
FFLinearCoefficients::usage = "FFLinearCoefficients[c1 x1 + c2 x2 + ..., {x1,x2,...}] returns the list of linear coefficients {c1,c2,...}.  Constant or higher-degree terms in the variables {x1,x2,...} are neglected."
FFSparseRowRules::usage = "FFSparseRowRules[c1 x1 + c2 x2 + ..., {x1,x2,...}] returns a list of rules of the form i -> ci,  where i is an integer and ci is the coefficient multiplying the i-th variable in {x1,x2,...}, for all non-vanishing coefficients.  Constant or higher-degree terms in the variables {x1,x2,...} are neglected."

FFDenseSolve::usage = "FFDenseSolve[eqs,vars] reconstructs and returns the solutions of the linear system eqs in the variables vars using a dense linear solver."
FFSparseSolve::usage = "FFSparseSolve[eqs,vars] reconstructs and returns the solutions of the linear system eqs in the variables vars using a sparse linear solver."
FFLinearFit::usage = "FFLinearFit[params, ansatz, rhs, vars, coeffs] reconstructs the coefficients coeffs, as functions of the free parameters params, which solve the linear fit equation ansatz == rhs.  The ansatz must be a linear expression of the form c1 f1 + c2 f2 + ... + f0, where {c1,c2,...} is the input list coeffs and {f1,f2,...} are rational functions of params and vars.  The ansatz can also be specified as a list {f1,f2,f3,...,f0} where f0 can be omitted if vanishing.  The input rhs is also a rational function of params and vars or a list of rational functions to be summed.  Both the functions f1,f2,... and rhs may also depend on additional auxiliary variables, which are rational functions of params and vars, to be defined via the option \"Substitutions\"."
FFInverse::usage = "FFInverse[mat] reconstructs the inverse of the square matrix mat via Gauss-Jordan elimination.  By default it uses a dense solver for the inversion.  The option \"Sparse\"->True can be passed to use a sparse solver instead."
FFTogether::usage = "FFTogether[expr], where expr is a rational expression, reconstructs a collected and GCD-simplified form of expr."

FFFunDeg::usage = "FFFunDeg[numdeg,dendeg] represents a generic rational functions with the total degree numdeg and dendeg for numerator and denominator respectively."

FFNewGraph::usage = "FFNewGraph[graphname] defines a new graph with name graphname.
FFNewGraph[graphname,inputnode,vars] is equivalent to FFNewGraph[graphname] followed by FFGraphInputVars[graphname,inputnode,vars]."
FFNewDummyGraph::usage = "FFNewDummyGraph[graphname,nparsin,nparsout] defines a dummy graph with nparsin input variables and nparsout output elements.  This graph cannot be evaluated but it can be used for functional reconstruction from stored evalutions."
FFDeleteGraph::usage = "FFDeleteGraph[graphname] deletes the graph with name graphname."
FFDeleteNode::usage = "FFDeleteNode[graph,node] deletes the specified node from a graph."
FFGraphInputVars::usage = "FFGraphInputVars[graph,nodename,vars] creates an input node representing the list of input variables vars.  Only the length of vars is relevant, and the names of the variables are not stored."
FFGraphOutput::usage = "FFGraphOutput[graph,node] selects the output node for a graph."
FFGraphDraw::usage = "FFGraphDraw[graph] returns a graphical representation of the specified dataflow graph.  With the option \"Pruned\"->True, only the nodes which are needed for evaluating the graph are drawn."
FFGraphEdges::usage = "FFGraphEdges[graph] returns a list of edges for the graph."

FFAlgDenseSolver::usage = "FFAlgDenseSolver[graph,node,{input},params,eqs,vars] defines a dense linear solver returning the solution of the equations eqs in the variables vars as function of the free input parameters params.  Only the solution for the non-vanishing variables is returned."
FFAlgSparseSolver::usage = "FFAlgSparseSolver[graph,node,{input},params,eqs,vars] defines a sparse linear solver returning the solution of the equations eqs in the variables vars as function of the free input parameters params."
FFAlgNodeDenseSolver::usage = "FFAlgNodeDenseSolver[graph,node,{input},neqs,vars] defines a dense linear solver returning the solution of the equations A[[i,1]] vars[[1]] + A[[i,2]] vars[[2]] + ... == b[[i]] in the variables vars, for i=1,...,neqs, where the entries returned by the input node are interpreted as the elements of the matrix (A b) in row-major order."
FFAlgSubgraphFit::usage = "FFAlgSubgraphFit[graph,node,{input},subgraph,vars,coeffs] is a subgraph algorithm which returns the solution for the coefficients coeffs, as functions of the input parameters, solving the linear fit problem coeffs[[1]] f[[1]] + coeffs[[2]] f[[2]] + ... + coeffs[[n-1]] f[[n-1]] == f[[n]], where f is the output list of subgraph, evaluated as a function of the list of variables vars concatenated with the input parameters represented by the input node.
FFAlgSubgraphFit[graph,node,{},subgraph,vars,coeffs] is a subgraph algorithm which returns the numerical solution for the coefficients coeffs, solving the linear fit problem coeffs[[1]] f[[1]] + coeffs[[2]] f[[2]] + ... + coeffs[[n-1]] f[[n-1]] == f[[n]], where f is the output list of subgraph, evaluated as a function of the variables vars."
FFAlgSubgraphReconstruct::usage = "FFAlgSubgraphReconstruct[graph,node,{input},subgraph,vars] reconstructs the output of subgraph as a list of functions of the variables vars and returns the coefficients of their monomials as functions of the input parameters.  The input of subgraph is the list of variables vars concatenated with the input parameters represented by the input node."
FFAlgSubgraphMultiFit::usage = "FFAlgSubgraphMultiFit[graph,node,inputs,subgraph,vars,take] is a subgraph node which performs several linear fits of the form c1 f1 + c2 f2 + ... == f0, where the elements {f1,f2,...,f0} are taken from the output of subgraph based on the pattern take."
FFSerializeSparseEqs::usage="FFSerializeSparseEqs[filename, params, eqs, vars, pattern, position] serializes the equations eqs in \"MX\" format in the file filename, formatted for communication with the finiteflow library."
FFAlgSerializedSparseSolver::usage="FFAlgSerializedSparseSolver[graph,node,{input},params,files,vars] defines a sparse solver using the equations serialized in files using FFSerializeSparseEqs, with respect to the variables vars, as a function of the free parameters params."
FFAlgDebug::usage = "FFAlgDebug[] shows a list of graphs and nodes defined in Mathematica, with their corresponding integer IDs, for debugging purposes."
FFAllAlgs::usage = "FFAllAlgs[] returns a list with all the algorithms in all the current graphs."
(*FFClearAlgs::usage = ""*)
FFAlgQ::usage = "FFAlgQ[graphname, nodename] returns True if the specified algorithm exists."
FFSolverResetNeededVars::usage = "FFSolverResetNeededVars[graph, node, vars, neededvars] redefines the set of needed variables of a dense or sparse linear system."
FFSolverOnlyHomogeneous::usage =  "FFSolverOnlyHomogeneous[graph, node] makes a linear solver return only the homogeneous part of its solution, i.e. without including the constant terms in the output."
FFSolverSparseOutput::usage = "FFSolverSparseOutput[graph, node] makes a sparse linear solver return a sparse representation of the solution matrix."
FFLearn::usage = "FFLearn[graph], executes the learning phase on the output node of graph."
FFSetLearningOptions::usage = "FFLearn[graph,node,options...] sets the learning options of the specified node in the graph."
FFLaurentLearn::usage = "FFLaurentLearn[graph] executes the learning phase on a Laurent expansion node, which must be the output node of graph.  It returns a list of two lists.  The first contains the starting power of the Laurent expansion of each element.  The second contains the order of the expansion requested for each element."
FFDenseSolverLearn::usage = "FFDenseSolverLearn[graph,vars] executes the learning phase on a dense solver or a linear fit, with unknowns vars, which must be the output node of graph."
FFSparseSolverLearn::usage = "FFSparseSolverLearn[graph,vars] executes the learning phase on a sparse solver, with unknowns vars, which must be the output node of graph."
FFDenseSolverGetInfo::usage = "FFDenseSolverGetInfo[graph,node,vars] returns the info about a dense solver or a linear fit, with unknowns vars, obtained during the learning step, namely dependent variables, independent variables and zero variables."
FFSparseSolverGetInfo::usage = "FFSparseSolverGetInfo[graph,node,vars] returns the info on a sparse solver, with unknowns vars, obtained during the learning step, namely dependent variables, independent variables and zero variables."
FFNonZeroesLearn::usage = "FFNonZeroesLearn[graph] executes the learning phase on a NonZero algorithm, which must be the output node of graph."
FFNonZeroesGetInfo::usage = "FFNonZeroesGetInfo[graph,ndode] returns the info on a NonZero algorithm obtained during the learning step."
FFMultiFitLearn::usage = "FFMultiFitLearn[graph,coefficients] executes the learning phase on a multi-fit subgraph algorithm, which must be the output node of graph."
FFMultiFitGetInfo::usage = "FFMultiFitGetInfo[graph,node,coefficients] returns the info on a multi-fit subgraph algorithm obtained during the learning step."
FFMultiFitSol::usage = "FFMultiFitSol[expr,learninfo], where expr represents the (reconstructed or symbolic) output of a multi-fit algorithm, and learninfo is the information obtained during its learning phase, formats expr as a list of lists of substitution rules representing the solutions of the fits."
FFSubgraphReconstructLearn::usage = "FFSubgraphReconstructLearn[graph,vars] executes the learning phase of a SubgraphReconstruct algorithm and returns, for each output element of the subgraph, a list of the form {{n1,n2,...},{d1,d2,...}} where {n1,n2,...} ({d1,d2,...}) is the list of non-vanishing monomials of its numerator (denominator) in the variables vars."
FFSubgraphReconstructGetInfo::usage = "FFSubgraphReconstructGetInfo[graph,node,vars] returns the info on a SubgraphReconstruct algorithm obtained during the learning step."
FFSubgraphReconstructSol::usage = "FFSubgraphReconstructSol[expr,learninfo], where expr represents the (reconstructed or symbolic) output of a SubgraphReconstruct algorithm and learninfo is the information obtained during its learning phase, returns a list of analytic expressions for the output of the subgraph."
FFAlgSimpleSubgraph::usage = "FFAlgSimpleSubgraph[graph,node,{input},subgraph] creates a simple subgraph node."
FFAlgMemoizedSubgraph::usage = "FFAlgMemoizedSubgraph[graph,node,{input},subgraph] creates a subgraph node which remembers the input and output of its last evaluation, avoiding repeated evaluations with the same input."
FFAlgSubgraphMap::usage = "FFAlgSubgraphMap[graph,node,inputs,subgraph] executes subgraph several times, using each list in inputs as input variables, and chains the outputs together."
FFSolverNIndepEqs::usage = "FFSolverNIndepEqs[graph,node] returns the number of independent equations of a linear system."
FFSolverIndepEqs::usage = "FFSolverIndepEqs[graph,node] returns a list of integers representing the indpendent equations of a linear system."
FFSparseSolverMarkAndSweepEqs::usage = "FFSparseSolverMarkAndSweepEqs[graph,node] executes the mark-and-sweep algorithm on a sparse system filtering out a subset of equations sufficient for returning the solution for the needed variables.   It returns the number of needed equations after the filtering."
FFSparseSolverDeleteUnneededEqs::usage = "FFSparseSolverDeleteUnneededEqs[graph,node] frees some memory by deleting the unneeded equations in a linear solver."
FFIndependentOf::usage = "FFIndependentOf[graph,varslist,var] returns True if the output of graph, as a function of the variables in varslist, is independent of the variable var, and False otherwise."
FFAlgRatFunEval::usage = "FFAlgRatFunEval[graph,node,{input},vars,funs] evaluates and returns the list of rational functions funs in the variables vars."
FFAlgRatExprEval::usage = "FFAlgRatExprEval[graph,node,{input},vars,expr] evaluates and returns the list of rational expressions expr in the variables vars.  The expressions do not need to be collected and will not be analytically expanded before numerical evaluation."
FFAlgRatFunEvalFromCoeffs::usage = "FFAlgRatFunEvalFromCoeffs[graph,node,{coeffinput,varsinput},coeffs,vars,funs] evaluates and returns the list of rational functions funs in the variables vars, where the coefficients of the monomials in the numerator and the denominator of the function are listed in coeffs and taken from the first input node coeffinput."
FFAlgRatNumEval::usage = "FFAlgRatFunEval[graph,node,ratnums] evaluates and returns the list of rational numbers ratnums."
FFAlgChain::usage = "FFAlgChain[graph,node,inputs] chains together the lists returned by its inputs."
FFAlgTake::usage="FFAlgTake[graph,node,inputs,takepattern] takes and returns selected elements from its inputs, as specified in takepattern."
FFAlgSlice::usage="FFAlgSlice[graph,node,{input},start,end], with integers start and end, returns the elements from position start to position end of its input.
FFAlgSlice[graph,node,{input},start], with integer start, returns the elements of its input starting from position start until the last one."
FFAlgAdd::usage = "FFAlgAdd[graph,node,inputs] adds the lists returned by the inputs element-wise and returns the result."
FFAlgMul::usage = "FFAlgMul[graph,node,inputs] multiplies the lists returned by the inputs element-wise and returns the result."
FFAlgMatMul::usage = "FFAlgMatMul[graph,node,{input1,input2},r1,c1,c2], with integers r1,c1,c2, interprets input1 and input2 as the elements of a r1 \[Times] c1 matrix and a c1 \[Times] c2 matrix respectively, in row-major order, and returns the result of the matrix multiplication input1.input2."
FFAlgSparseMatMul::usage = "FFAlgSparseMatMul[graph,node,{input1,input2},r1,c1,c2,nonzerocols1,nonzerocols2] is analogous to FFAlgMatMul[graph,node,{input1,input2},r1,c1,c2] except that input1 and input2 only return the potentially non-vanishing matrix elements of the inputs.  The arguments nonzerocols1 and nonzerocols2 are lists of lists with the potentially non-vanishing columns in each row for the two input matrices respectively."
FFAlgNonZeroes::usage = "FFAlgNonZeroes[graph,node,{input}] returns the non-vanishing elements of its input."
FFAlgTakeAndAdd::usage="FFAlgTakeAndAdd[graph,node,inputs,takeelements] is a TakeAndAdd algorithm which takes lists of selected elements from its inputs, according to the takeelements pattern, and adds all elements in each list."
FFTotalDegrees::usage="FFTotalDegrees[graph] computes and returns the total degrees of each entry of the output of a graph."
FFVarsDegrees::usage="FFVarsDegrees[graph] computes and internally stores the degrees with respect to each variable, for all the outputs of a graph."
FFAllDegrees::usage="FFAllDegrees[graph] computes the total degrees, as well as the partial degrees with respect to each variable, for all the outputs of a graph.  The total degrees are also returned."
FFSample::usage="FFSample[graph] evaluates a graph for a set of sample points, which depends on the reconstruction options."
FFReconstructFromCurrentEvaluations::usage="FFReconstructFromCurrentEvaluations[graph,vars] attempts to analytically reconstruct the output of a graph using the numerical evaluations which have already been performed and stored."
FFReconstructNumeric::usage="FFReconstructNumeric[graph] performs a numerical reconstruction over the rational field of the output of a graph with no input node."
FFReconstructFromCurrentEvaluationsMod::usage="FFReconstructFromCurrentEvaluationsMod[graph,vars] is the same as FFReconstructFromCurrentEvaluations[graph,vars] but the reconstruction is performed modulo the prime specified in the options."
FFMissingPoints::usage = "Returned when there are not enough sample points for reconstructing a function on a given prime field."
FFMissingPrimes::usage = "Returned when sample points from additional prime fields are needed for reconstructing a function."

FFRegisterAlgorithm::usage = "Low level interface for implementing new native algorithms."

FFSparseEqsToJSON::usage="FFSparseEqsToJSON[outputfilename,params,eqs,vars,pattern,position] serializes the equations eqs in the variables vars, and depending of the free parameters params, in JSON format."
FFSparseSystemToJSON::usage="FFSparseSystemToJSON[outputfilename,neqs,vars,pars,filelist] creates a JSON file with the information about a sparse system of equations in the variables vars, and depending of the free parameters params, which have been serialized in JSON format in the files of the list filelist.  These JSON files can be created using FFSparseEqsToJSON."
FFAlgJSONSparseSolver::usage = "FFAlgJSONSparseSolver[graph,node,{input},filename] defines a sparse solver for a system defined in the JSON file filename.  The JSON file can be created using FFSparseSystemToJSON."

FFRatFunToJSON::usage="FFRatFunToJSON[outputfilename,vars,funcs] exports the list of rational functions funcs in the variables vars in JSON format."
FFAlgJSONRatFunEval::usage="FFAlgJSONRatFunEval[graph,node,{input},filename] creates a node evaluating the list of rational functions serialized in JSON format in the file filename.  The JSON file can be created using FFRatFunToJSON."
FFRatFunFromJSON::usage="FFRatFunFromJSON[filename,vars] unserialized and returns a list of rational functions in the variables vars, stored in JSON format in the file filename."

FFListSubsetToJSON::usage = "FFListSubsetToJSON[outputfilename,list,subset] defines a subset of the input list and serialized it JSON format."
FFListSubsetFromJSON::usage "FFListSubsetFromJSON[filename,list] returns a subset of list, specified by the JSON file filename."

FFAlgLinearFit::usage = "FFAlgLinearFit[graph,node,{input},params,ansatz, rhs, vars, coeffs] returns a node performing a linear fit.  See the usage FFLinearFit for more information."

FFAlgLaurent::usage = "FFAlgLaurent[graph,node,{input},subgraph,order] defines a subgraph node which computes the coefficients of the Laurent expansion of the output of subgraph around x=0, where x is the first input variable of subgraph.  If order is an integer, the expansion is truncated at x^order for all output elements.  If order is a list, the expansion of the i-th output element is truncated at x^order[[i]]."
FFLaurentSol::usage = "FFLaurentSol[expr,x,learninfo] formats the output of a Laurent expansion node as series expansion in x, using the information learninfo returned by FFLearn."

FFReconstructFunction::usage = "FFReconstructFunction[graph,vars] analytically reconstructs the output of graph as a rational function in the variables vars, by performing several numerical evaluations and reconstructing analytic expressions from these.  Note that for univariate functions, multi-threading is not used during the evaluation of the graph (for that, FFParallelReconstructUnivariate can be used)."
FFParallelReconstructUnivariate::usage = "FFParallelReconstructUnivariate[graph,{x}] reconstructs the output of graph as a univariate rational function in the variable x.  Numerical evaluations are performed in parallel, and additional sample points are added until the reconstruction is successful."
FFReconstructFunctionMod::usage = "FFReconstructFunctionMod[graph,vars] is the same as FFReconstructFunction[graph,vars] but reconstruction is performed modulo the prime specified in the options."

FFDenseSolverSol::usage = "FFDenseSolverSol[expr,learninfo], where expr represents the (reconstructed or symbolic) output of a dense solver or a linear fit, and learinfo is the information obtained during its learning phase, formats expr as list of substitution rules representing the solution of the system."
FFSparseSolverSol::usage = "FFSparseSolverSol[expr,learninfo], where expr represents the (reconstructed or symbolic) output of a sparse solver, and learinfo is the information obtained during its learning phase, formats expr as list of substitution rules representing the solution of the system."
FFNonZeroesSol::usage = "FFNonZeroesSol[expr,learninfo], where the list expr is the output of a NonZero algorithm and learninfo is the information returned by its learning phase, puts back the zero elements in the list expr."

FFNThreads::usage = "Default number of threads used (default: Automatic)"

FFAutomaticNThreads::usage = "FFAutomaticNThreads[] returns the default number of threads used when Automatic is selected"

FFNParsOut::usage = "FFNParsOut[graph] returns the length of the output of a graph.
FFNParsOut[graph,node] returns the length of the output list of a node in a graph."

FFMakeMutable::usage = "FFMakeMutable[graph,node] makes a node mutable, if possible."
FFGraphPrune::usage = "FFGraphPrune[graph] deletes all the nodes of graph which are not needed for the evaluation of its output node."

FFDumpDegrees::usage="FFDumpDegrees[graph,outputfilename] serializes the information on the total and partial degrees of the output of graph.  It must be called after FFAllDegrees."
FFLoadDegrees::usage="FFDumpDegrees[graph,filename] loads the information on total and partial degrees of graph, as stored in the file filename."
FFDumpSamplePoints::usage="FFDumpSamplePoints[graph,outputfilename] serializes a list of sample points, i.e. a list of inputs at which graph needs to be evaluated in order to reconstruct the analytic expression of its output."
FFDumpEvaluations::usage="FFDumpEvaluations[graph,outputfilename] serializes the currently stored evaluations of graph."
FFLoadEvaluations::usage="FFLoadEvaluations[graph,filelist] loads the evaluations of graph serialized in the files listed in filelist."
FFSamplesFileSize::usage="FFSamplesFileSize[filename] returns the number of sample points listed in the file filename."
FFNParsFromDegreesFile::usage="FFNParsFromDegreesFile[filename] returns the number of input parameters of a graph, reading it from the file filename which stores the information on its degrees."
FFSampleFromPoints::usage="\
FFSampleFromPoints[graph,filename] evaluates the graph at the points stored in filename.  The evaluations are parallelized over an automatically chosen number of threads.
FFSampleFromPoints[graph,filename,nthreads] is equivalent to FFSampleFromPoints[graph,filename], but the evaluations are parallelized over nthreads threads.
FFSampleFromPoints[graph,filename,start,npoints] evaluates the graph at a contiguous subset of npoints sample points, taken from the ones stored in filename, starting from the one at position start (counting from zero).  The evaluations are parallelized over an automatically chosen number of threads.
FFSampleFromPoints[graph,filename,start,npoints,nthreads] is equivalent to FFSampleFromPoints[graph,filename,start,npoints], but the evaluations are parallelized over nthreads threads."
FFNSamplePoints::usage = "FFNSamplePoints[graph] returns a list of length two.  The first element is the total number of sample points needed for recostructing the full output of graph.  The second element is a list of integers representing the number of sample points needed for the reconstruction of each element of the output of graph."

FFGraphEvaluate::usage="FFGraphEvaluate[graph,point] evaluates graph at point, where point is a list of integers.  The prime field may be changed passing the option \"PrimeNo\"."
FFGraphEvaluateMany::usage="FFGraphEvaluateMany[graph,points] evaluates graph at the specified list of points.  The prime field may be changed globally using the option \"PrimeNo\", or individually for each point by appending an additional entry with the index of the prime to be used.  By default, evaluations are performed in parallel."
FFPrimeNo::usage="FFPrimeNo[i] with i>=0 returns the i-th hardcoded prime used by finiteflow."
FFMulInv::usage="FFMulInv[z,p] returns the multiplicative inverse of the integer z module a prime p."
FFRatMod::usage="FFRatMod[z,p] returns z mod p, where z is a rational number and p is a prime."

FFRatRec::usage="FFRatRec[a,n], where a and n are integers, returns a rational q such that q mod n = a, computed using Wang's rational reconstruction algorithm.
FFRatRec[{a1,a2,...},n] is equivalent to {FFRatRec[a1],n],FFRatRec[a2],n],...}."


FF::badrational = "Argument `1` is not a rational number."
FF::badfun = "Argument is not a polynomial or a rational function in the specified variables `1` with rational coefficients."
FF::badfunarg = "Argument `1` is not a polynomial or a rational function in the specified variables `2` with rational coefficients."
FF::badfuncoeff = "Argument `1` is not a polynomial or a rational function in the specified variables `2` with the specified coefficients."
FF::badvars = "`1` is not a non-empty list of variables."
FF::badsystem = "Argument is not a list of equalities."
FF::badbooleanflag = "The option `1` must be True, False or Automatic."
FF::baduintflag = "`1` must be a non-negative integer or Automatic."
FF::nolinear = "The system is not linear in the variables `1`."
FF::noint = "`1` is not an integer."
FF::noint32 = "`1` is not an integer in the 32-bit range."
FF::noint32list = "`1` is not a list of integers in the 32-bit range."
FF::nouint32list = "`1` is not a list of positive integers in the 32-bit range."
FF::noint64 = "`1` is not an integer in the 64-bit range."
FF::noint64list = "`1` is not a list of integers in the 64-bit range."
FF::nolib = "fflow library cannot be loaded: try setting your $LibraryPath, and then call FFLoadLib[]."
FF::badsquaremat = "The input is not a square matrix."
FF::badrule = "`1` is not a list of substitution rules."
FF::badvarule = "Invalid variable on l.h.s. of substitution rule."
FF::noregfun = "No registered function with identifier `1`."
FF::badregfunvars = "The registered expression `1` depends on `2` variables, but `3` are required."
FF::badpattern = "Variables `1` match the variables pattern but they are not in the list of unknowns."

FF::noalg = "No algorithm with identifier `1`."
FF::nograph = "No graph with identifier `1`."
FF::badalgvars = "The algorithm `1` depends on `2` variables, but `3` are required."
FF::badneededvars = "Needed variables should be a subset of the unknowns."

FF::nonratsub = "Found invalid subexpression `1`"


Begin["`Private`"]


FFMulInv[a_,p_]:=Mod[ExtendedGCD[a,p][[2,1]],p];
FFRatMod[a_,p_]:=Mod[Numerator[a]*FFMulInv[Denominator[a],p],p];


(* Constants *)
Unprotect[FFInt64Max];

FFInt32Max = 2^31-1;
FFInt64Max = 2^63-1;

Protect[FFInt64Max];


FFRationalQ[x_Rational] := True;
FFRationalQ[x_Integer] := True;
FFRationalQ[x_] := False;


CheckedInt[a_] := If[IntegerQ[a], a, Message[FF::noint, a]; Throw[$Failed]];
CheckedInt32[a_] := If[IntegerQ[a] && (Abs[a]<=FFInt32Max), a, Message[FF::noint, a]; Throw[$Failed]];
CheckedInt32Range[a_] := If[Abs[a]<=FFInt32Max, a, Message[FF::noint32, a]; Throw[$Failed]];
CheckedInt64[a_] := If[IntegerQ[a] && (Abs[a]<=FFInt64Max), a, Message[FF::noint, a]; Throw[$Failed]];
CheckedInt64Range[a_] := If[Abs[a]<=FFInt64Max, a, Message[FF::noint64, a]; Throw[$Failed]];

Int32ListError[a_] := (Message[FF::noint32list, a]; Throw[$Failed]);
CheckedInt32List[a_] := Int32ListError[a];
CheckedInt32List[a_List] := If[TrueQ[And@@((IntegerQ[#] && Abs[#]<=FFInt32Max)&/@a)], a, Int32ListError[a]];

UInt32ListError[a_] := (Message[FF::nouint32list, a]; Throw[$Failed]);
CheckedUInt32List[a_] := UInt32ListError[a];
CheckedUInt32List[a_List] := If[TrueQ[And@@((IntegerQ[#] && #>=0 && #<=FFInt32Max)&/@a)], a, UInt32ListError[a]];

Int64ListError[a_] := (Message[FF::noint64list, a]; Throw[$Failed]);
CheckedInt64List[a_] := Int64ListError[a];
CheckedInt64List[a_List] := If[TrueQ[And@@((IntegerQ[#] && Abs[#]<=FFInt64Max)&/@a)], a, Int64ListError[a]];

CheckVariables[vars_]:=If[!TrueQ[(vars[[0]] == List) && (Length[vars] > 0)],
                           Message[FF::badfun, vars]; Throw[$Failed]
	                   ];


FFCoefficientRules[expr_,vars_] := (((#[[1]]-1)->#[[2]])&/@(ArrayRules[SparseArray[CoefficientList[expr,vars]]][[;;-2]]));
FFCoefficientRules[0,vars_] := {};


FFLinearCoefficients[0,vars_]:=ConstantArray[0,Length[vars]];
FFLinearCoefficients[expr_,vars_]:=Normal[CoefficientArrays[expr,vars][[2]]];
FFSparseRowRules[0,columnvec_]:={};
FFSparseRowRules[rowexpr_,columnvec_]:=#[[1,1]]->#[[2]]&/@SortBy[ArrayRules[CoefficientArrays[rowexpr,columnvec][[2]]][[;;-2]],First];


PolyCoefficientRules[poly_,vars_] :=  If[AllTrue[#,FFRationalQ[#[[2]]]&],
                                         #,
                                         Message[FF::badfunarg, poly, vars]; Throw[$Failed]
                                       ]&[FFCoefficientRules[poly,vars]];


PolyCoefficientRulesCoeffMap[poly_,vars_,map_] :=  If[SubsetQ[Keys[map],(#[[2]]&/@#)],
                                                    #,
                                                    Message[FF::badfuncoeff, poly, vars]; Throw[$Failed]
                                                    ]&[FFCoefficientRules[poly,vars]];


LinearEqCoeffs[expr_, vars_, applyfun_] := Module[
    {res, nvars},
    nvars = Length[vars];
    res = applyfun[Normal@CoefficientArrays[expr, vars]];
    If[TrueQ[Length[res]>2], Message[FF::nolinear, vars]; Throw[$Failed]];
    If[TrueQ[Length[res]==1],
      Join[ConstantArray[0,nvars],{-res[[1]]}],
      Join[res[[2]],{-res[[1]]}]
    ]
];


SparseLinearEqCoeffs[expr_, vars_, applyfun_] := Module[
    {res, nvars, ret},
    nvars = Length[vars];
    res = applyfun[(ArrayRules[#][[;;-2]])&/@(CoefficientArrays[{expr}, vars])];
    If[TrueQ[Length[res]>2], Message[FF::nolinear, vars]; Throw[$Failed]];
    If[TrueQ[Length[res]==1],
      If[Length[res[[1]]]==0,{{},{}},{{Length[vars]},{-res[[1,1,2]]}}],
      res[[2]] = SortBy[res[[2]], First];
      ret = {(#[[1,2]]-1)&/@res[[2]],#[[2]]&/@res[[2]]};
      If[Length[res[[1]]]!=0, ret = {Join[ret[[1]],{Length[vars]}], Join[ret[[2]],{-res[[1,1,2]]}]}];
      ret
    ]
];


SparseLinearEqCoeffsWithPattern[expr_, totnvars_, pattern_, position_, applyfun_] := Module[
    {res, ret, vars},
    vars = pattern[expr];
    res = applyfun[(ArrayRules[#][[;;-2]])&/@(CoefficientArrays[{expr}, vars])];
    If[TrueQ[Length[res]>2], Message[FF::nolinear, vars]; Throw[$Failed]];
    If[TrueQ[Length[res]==1],
      If[Length[res[[1]]]==0,{{},{}},{{totnvars},{-res[[1,1,2]]}}],
      res[[2]] = SortBy[res[[2]], position[vars[[#[[1,2]]]]]&];
      ret = {position[vars[[#[[1,2]]]]]&/@res[[2]],#[[2]]&/@res[[2]]};
      If[Length[res[[1]]]!=0, ret = {Join[ret[[1]],{totnvars}], Join[ret[[2]],{-res[[1,1,2]]}]}];
      ret
    ]
];


toFFInternalBooleanFlag[var_, Automatic]:=0;
toFFInternalBooleanFlag[var_, True]:=1;
toFFInternalBooleanFlag[var_, False]:=-1;
toFFInternalBooleanFlag[var_, other_]:=(Message[FF::badbooleanflag, var]; Throw[$Failed]);


toFFInternalUnsignedFlag[var_, Automatic]:=-1;
toFFInternalUnsignedFlag[var_, n_Integer]:=If[TrueQ[n>=0], n, Message[FF::baduintflag, var]; Throw[$Failed]];
toFFInternalUnsignedFlag[var_, other_]:=(Message[FF::baduintflag, var]; Throw[$Failed]);


toFFInternalPoly[poly_,vars_] := ({#[[1]],ToString[#[[2]],InputForm]})&/@PolyCoefficientRules[poly, vars];
toFFInternalRatFun[ratfun_,vars_] := {toFFInternalPoly[Numerator[ratfun],vars],toFFInternalPoly[Denominator[ratfun],vars]};
toFFInternalRatFun[0,vars_] := {Length[vars]}; (* optimization \[Rule] this represents a vanishing function *)


toFFInternalPolyCoeffMap[poly_,vars_,coeffmap_] := ({#[[1]],coeffmap[#[[2]]]})&/@PolyCoefficientRulesCoeffMap[poly, vars, coeffmap];
toFFInternalRatFunCoeffMap[ratfun_,vars_,coeffmap_] := {toFFInternalPolyCoeffMap[Numerator[ratfun],vars,coeffmap],toFFInternalPolyCoeffMap[Denominator[ratfun],vars,coeffmap]};
toFFInternalRatFunCoeffMap[0,vars_,coeffmap_] := {Length[vars]}; (* optimization \[Rule] this represents a vanishing function *)


fromFFInternalPoly[ipoly_,vars_] := Plus@@((ToExpression[#[[2]]]Times@@(vars^#[[1]]))&/@ipoly);
fromFFInternalRatFun[ifun_,vars_] := fromFFInternalPoly[ifun[[1]],vars]/fromFFInternalPoly[ifun[[2]],vars];

fromFFInternalPoly[$Failed,vars_] := $Failed;
fromFFInternalRatFun[$Failed,vars_] := $Failed;
fromFFInternalRatFun[FFBadShift,vars_] := FFBadShift;

fromFFInternalRatFun[ifun_FFFunDeg,vars_] := ifun;
fromFFInternalRatFun[0,vars_] := 0;


toFFInternalPolyPoly[poly_,tauvars_,params_] := ({#[[1]],toFFInternalPoly[#[[2]],params]})&/@FFCoefficientRules[poly, tauvars];
toFFInternalPolyRatFun[ratfun_,tauvars_,params_] := {toFFInternalPolyPoly[Numerator[ratfun],tauvars,params],toFFInternalPolyPoly[Denominator[ratfun],tauvars,params]};


toFFJSONPoly[poly_,vars_] := ({ToString[#[[2]],InputForm],#[[1]]})&/@PolyCoefficientRules[poly, vars];
toFFJSONRatFun[ratfun_,vars_] := {Length[#],#}&/@{toFFJSONPoly[Numerator[ratfun],vars],toFFJSONPoly[Denominator[ratfun],vars]};
toFFJSONRatFun[0,vars_] := 0;


FFAutomaticNThreads[]:=FFDefaultNThreadsImplem[];


(*FFMulInv[a_Integer, mod_Integer] := Catch[FFMulInvImplem[CheckedInt64Range[a], CheckedInt64Range[mod]]];*)


If[!TrueQ[FFAlreadyLoaded],
  FFGraphId = Association[{}];
  FFAlgId = Association[{}];
  FFGraphInputs = Association[{}];
  FFNThreads = Automatic;
];

FFAlgDebug[]:=Print[{FFGraphId, FFAlgId}];

FFAllAlgs[]:=Keys[FFAlgId];

(*FFClearAlgs[]:=(FFClearAllAlgsImplem[]; FFAlgId = Association[{}]; FFAlgVars = Association[{}];);*)

FFGraphQ[graphid_] := IntegerQ[FFGraphId[graphid]];
FFAlgQ[graphid_,id_] := IntegerQ[FFAlgId[{graphid,id}]];

GetGraphId[gid_]:=(If[!FFGraphQ[gid], Message[FF::nograph,gid]; Throw[$Failed];];
                   FFGraphId[gid]);
GetAlgId[gid_, id_]:=(If[!FFAlgQ[gid,id], Message[FF::noalg,gid->id]; Throw[$Failed];];
                      FFAlgId[{gid,id}]);

FFValidateId[graphid_,algid_]:=FFGrapQ[graphid];
FFValidateId[any___]:=False;

FFRemoveGraphKeys[graphid_]:=(KeyDropFrom[FFGraphId,graphid];
                                 KeyDropFrom[FFAlgId,Select[FFAllAlgs[],TrueQ[#[[1]]==graphid]&]];);
FFDeleteGraph[graphid_]:=If[FFGraphQ[graphid], FFDeleteGraphImplem[FFGraphId[graphid]]; FFRemoveGraphKeys[graphid];];
FFNewGraph[graphid_]:=(If[FFGraphQ[graphid], FFDeleteGraph[graphid]]; FFGraphId[graphid]=FFNewGraphImplem[]);
FFNewDummyGraph[graphid_,nparsin_,nparsout_]:=Catch[(If[FFGraphQ[graphid], FFDeleteGraph[graphid]]; FFGraphId[graphid]=FFNewDummyGraphImplem[CheckedInt32[nparsin],CheckedInt32[nparsout]])];
FFNewGraph[graphid_,vars_List]:=(FFNewGraph[graphid]; FFGraphInputVars[graphid,vars]);
FFNewGraph[graphid_,inputnode_,vars_List]:=(FFNewGraph[graphid]; FFGraphInputVars[graphid,inputnode,vars]);
FFDeleteNode[graphid_,nodeid_]:=Module[{ret},
  ret = $Failed;
  If[FFAlgQ[graphid,nodeid],
    ret = FFDeleteNodeImplem[FFGraphId[graphid],FFAlgId[{graphid,nodeid}]];
    If[!TrueQ[ret==$Failed], KeyDropFrom[FFAlgId,{{graphid,nodeid}}]];
  ];
  ret
];

FFGraphInputVars[graphid_,vars_List]:=Module[{ret},
  ret = FFSetInputVarsImplem[GetGraphId[graphid],Length[vars]];
  If[TrueQ[ret==$Failed], Return[$Failed]];
  KeyDropFrom[FFAlgId,{{graphid,FFGraphInputs[graphid]}}];
  KeyDropFrom[FFGraphInputs,graphid];
  FFAlgId[{graphid,vars}]=0;
  FFGraphInputs[graphid]=vars;
  True
];
FFGraphInputVars[graphid_,shortcut_,vars_List]:=Module[{ret},
  ret = FFSetInputVarsImplem[GetGraphId[graphid],Length[vars]];
  If[TrueQ[ret==$Failed], Return[$Failed]];
  KeyDropFrom[FFAlgId,{{graphid,FFGraphInputs[graphid]}}];
  KeyDropFrom[FFGraphInputs,graphid];
  FFAlgId[{graphid,shortcut}]=0;
  FFGraphInputs[graphid]=shortcut;
  True
];

FFRegisterAlgorithm[algregfun_, gid_, id_, inputs_, args_List]:=Module[
    {present,idno},
     If[!TrueQ[FFGraphQ[gid]], Message[FF::nograph,gid]; Return[$Failed]];
     If[!AllTrue[inputs, FFAlgQ[gid,#]&],Message[FF::noalg,gid->SelectFirst[inputs,!FFAlgQ[gid,#]&]]; Return[$Failed]];
     If[FFAlgQ[gid,id], If[FFDeleteNode[gid,id]==$Failed, Return[$Failed]]];
     idno = algregfun[FFGraphId[gid],FFAlgId[{gid,#}]&/@inputs,args];
     If[TrueQ[idno == $Failed], Return[$Failed]];
     (*FFRegisterImplem[toFFInternalRatFun[expr,vars]];*)
     FFAlgId[{gid,id}] = idno;
     True
];

FFGraphOutput[graphid_,nodeid_]:=FFGraphSetOutputImplem[GetGraphId[graphid],GetAlgId[graphid,nodeid]];


FFNParsOut[gid_]:=Catch[FFGraphNParsOutImplem[GetGraphId[gid]]];
FFNParsOut[gid_,id_]:=Catch[FFNodeNParsOutImplem[GetGraphId[gid],GetAlgId[gid,id]]];


FFMakeMutable[gid_,id_]:=Catch[FFAlgMakeNodeMutableImplem[GetGraphId[gid],GetAlgId[gid,id]]];


Options[FFGraphEdges]:={"Pruned"->False};
FFGraphEdges[graph_,OptionsPattern[]]:=Module[{ret,map},
  Catch[
    ret = FFGraphEdgesImplem[GetGraphId[graph],toFFInternalBooleanFlag["Pruned",OptionValue["Pruned"]]];
    If[TrueQ[ret == $Failed],Throw[$Failed]];
    If[TrueQ[ret == {}],Throw[{}]];
    
    map = Association[{}];
    (map[FFAlgId[#]]=#[[2]];)&/@Select[FFAllAlgs[],First[#]==graph&];
    ret = map/@ret;
    (Rule@@#)&/@ArrayReshape[ret,{Length[ret]/2,2}]
  ]
];
Options[FFGraphDraw]=Options[FFGraphEdges];
FFGraphDraw[graph_,opt:OptionsPattern[]]:=
  If[TrueQ[#==$Failed],
     $Failed,
     LayeredGraphPlot[#,
                   DirectedEdges -> True,
                   VertexLabeling -> True(*,
                   VertexRenderingFunction -> ({White, EdgeForm[Black], Disk[#, .33], Black, Text[#2, #1]} &)*)]
     ]&@FFGraphEdges[graph,opt];


Options[FFGraphNodes]:={"Pruned"->False};
FFGraphNodes[graph_,OptionsPattern[]]:=Module[{ret,map},
  Catch[
    ret = FFGraphNodesImplem[GetGraphId[graph],toFFInternalBooleanFlag["Pruned",OptionValue["Pruned"]]];
    If[TrueQ[ret == $Failed],Throw[$Failed]];
    If[TrueQ[ret == {}],Throw[{}]];
    
    map = Association[{}];
    (map[FFAlgId[#]]=#[[2]];)&/@Select[FFAllAlgs[],First[#]==graph&];
    map/@ret
  ]
];


FFGraphPrune[gid_]:=Module[
 {ret,nodes,newnodes},
  Catch[
    ret = FFGraphPruneImplem[GetGraphId[gid]];
    If[TrueQ[ret==$Failed],Throw[$Failed]];
    nodes = Select[FFAllAlgs[],First[#]==gid&];
    newnodes = FFGraphNodes[gid];
    If[TrueQ[newnodes==$Failed],Throw[$Failed]];
    If[!TrueQ[MemberQ[newnodes,#[[2]]]],
      KeyDropFrom[FFAlgId,{#}]
    ]&/@nodes;
  ]
];


FFReconstructOptions = {"Checks", "UChecks", "MaxSingularPoints", "StartingPrimeNo", "MaxPrimes", "MaxDegree", "PrintDebugInfo"};
Options[FFAlgorithmSetReconstructionOptions]=Options[FFAlgorithmSetDefaultReconstructionOptions]=(#->Automatic)&/@FFReconstructOptions;
FFAlgorithmSetReconstructionOptions[OptionsPattern[]]:=toFFInternalUnsignedFlag[#,OptionValue[#]]&/@FFReconstructOptions;


Options[FFTotalDegrees]:=Options[FFAlgorithmSetReconstructionOptions];
FFTotalDegrees[gid_,opt:OptionsPattern[]]:=FFTotalDegreesImplem[GetGraphId[gid],FFAlgorithmSetReconstructionOptions[opt]];


Options[FFVarsDegrees]:=Options[FFAlgorithmSetReconstructionOptions];
FFVarsDegrees[gid_,opt:OptionsPattern[]]:=FFVarsDegreesImplem[GetGraphId[gid],FFAlgorithmSetReconstructionOptions[opt]];


Options[FFAllDegrees]:=Options[FFAlgorithmSetReconstructionOptions];
FFAllDegrees[gid_, nthreads_:FFNThreads, opt:OptionsPattern[]]:=FFAllDegreesImplem[GetGraphId[gid],toFFInternalUnsignedFlag["nthreads", nthreads], FFAlgorithmSetReconstructionOptions[opt]];
FFAllDegrees[gid_, opt:OptionsPattern[]]:=FFAllDegrees[gid,FFNThreads,opt];


Options[FFSample]:=Options[FFAlgorithmSetReconstructionOptions];
FFSample[gid_, nthreads_:FFNThreads, opt:OptionsPattern[]]:=FFSampleImplem[GetGraphId[gid],toFFInternalUnsignedFlag["nthreads", nthreads], FFAlgorithmSetReconstructionOptions[opt]];
FFSample[gid_, opt:OptionsPattern[]]:=FFSample[gid,FFNThreads,opt];


Options[FFSampleFromPoints]:=Options[FFAlgorithmSetReconstructionOptions];
FFSampleFromPoints[gid_, file_String, nthreads_:FFNThreads, opt:OptionsPattern[]]:=FFSampleFromPointsImplem[GetGraphId[gid],toFFInternalUnsignedFlag["nthreads", nthreads], FFAlgorithmSetReconstructionOptions[opt],
                                                                                                                                       file,0,0];
FFSampleFromPoints[gid_, file_String, start_, size_, nthreads_:FFNThreads, opt:OptionsPattern[]]:=Catch[FFSampleFromPointsImplem[GetGraphId[gid],toFFInternalUnsignedFlag["nthreads", nthreads], FFAlgorithmSetReconstructionOptions[opt],
                                                                                                                                                            file,CheckedInt32[start],CheckedInt32[size]]];
FFSampleFromPoints[gid_, file_String, opt:OptionsPattern[]]:=FFSampleFromPoints[gid,file,FFNThreads,opt];
FFSampleFromPoints[gid_, file_String, start_, size_, opt:OptionsPattern[]]:=FFSampleFromPoints[gid,file,start,size,FFNThreads,opt];


Options[FFNSamplePoints]:=Options[FFAlgorithmSetReconstructionOptions];
FFNSamplePoints[gid_, file_String, nthreads_:FFNThreads, opt:OptionsPattern[]]:=Catch[FFNSamplePointsImplem[GetGraphId[gid],toFFInternalUnsignedFlag["nthreads", nthreads], FFAlgorithmSetReconstructionOptions[opt], file]];
FFNSamplePoints[gid_, file_String, opt:OptionsPattern[]]:=FFNSamplePoints[gid,file,FFNThreads,opt];
FFNSamplePoints[gid_, nthreads_:FFNThreads, opt:OptionsPattern[]]:=Catch[FFNSamplePointsImplem[GetGraphId[gid],toFFInternalUnsignedFlag["nthreads", nthreads], FFAlgorithmSetReconstructionOptions[opt], ""]];
FFNSamplePoints[gid_, opt:OptionsPattern[]]:=FFNSamplePoints[gid,FFNThreads,opt];


Options[FFReconstructFromCurrentEvaluations]:=Options[FFAlgorithmSetReconstructionOptions];
FFReconstructFromCurrentEvaluations[gid_,vars_, nthreads_:FFNThreads, opt:OptionsPattern[]]:=Module[
  {res},
  res = FFReconstructFromCurrentEvaluationsImplem[GetGraphId[gid], toFFInternalUnsignedFlag["nthreads", nthreads], FFAlgorithmSetReconstructionOptions[opt]];
  If[!TrueQ[res[[0]]==List], Return[res]];
  fromFFInternalRatFun[#,vars]&/@res
];
FFReconstructFromCurrentEvaluations[gid_,vars_, opt:OptionsPattern[]]:=FFReconstructFromCurrentEvaluations[gid,vars,FFNThreads,opt];


Options[FFReconstructFromCurrentEvaluationsMod]:=Options[FFAlgorithmSetReconstructionOptions];
FFReconstructFromCurrentEvaluationsMod[gid_,vars_, nthreads_:FFNThreads, opt:OptionsPattern[]]:=Module[
  {res},
  res = FFReconstructFromCurrentEvaluationsModImplem[GetGraphId[gid], toFFInternalUnsignedFlag["nthreads", nthreads], FFAlgorithmSetReconstructionOptions[opt]];
  If[!TrueQ[res[[0]]==List], Return[res]];
  fromFFInternalRatFun[#,vars]&/@res
];
FFReconstructFromCurrentEvaluationsMod[gid_,vars_, opt:OptionsPattern[]]:=FFReconstructFromCurrentEvaluationsMod[gid,vars,FFNThreads,opt];


Options[FFReconstructUnivariate]:=Options[FFAlgorithmSetReconstructionOptions];
FFReconstructUnivariate[gid_,vars_, opt:OptionsPattern[]]:=Module[
  {res},
  res = FFReconstructUnivariateImplem[GetGraphId[gid], FFAlgorithmSetReconstructionOptions[opt]];
  If[!TrueQ[res[[0]]==List], Return[res]];
  fromFFInternalRatFun[#,vars]&/@res
];


Options[FFReconstructUnivariateMod]:=Options[FFAlgorithmSetReconstructionOptions];
FFReconstructUnivariateMod[gid_,vars_, opt:OptionsPattern[]]:=Module[
  {res},
  res = FFReconstructUnivariateModImplem[GetGraphId[gid], FFAlgorithmSetReconstructionOptions[opt]];
  If[!TrueQ[res[[0]]==List], Return[res]];
  fromFFInternalRatFun[#,vars]&/@res
];


FFToExpression[0]=0;
FFToExpression[a_]:=ToExpression[a];


Options[FFReconstructNumeric]:=Options[FFAlgorithmSetReconstructionOptions];
FFReconstructNumeric[gid_, opt:OptionsPattern[]]:=Module[
  {res},
  res = FFReconstructNumericImplem[GetGraphId[gid], FFAlgorithmSetReconstructionOptions[opt]];
  If[!TrueQ[res[[0]]==List], Return[res]];
  FFToExpression/@res
];


FFDumpDegrees[gid_,file_String]:=FFDumpDegreesImplem[GetGraphId[gid],file];
FFLoadDegrees[gid_,file_String]:=FFLoadDegreesImplem[GetGraphId[gid],file];
Options[FFDumpSamplePoints]:=Options[FFAlgorithmSetReconstructionOptions];
FFDumpSamplePoints[gid_,file_String,opt:OptionsPattern[]]:=FFDumpSamplePointsImplem[GetGraphId[gid],file, FFAlgorithmSetReconstructionOptions[opt]];
FFDumpEvaluations[gid_,file_String]:=FFDumpEvaluationsImplem[GetGraphId[gid],file];
FFLoadEvaluations[gid_,files_List]:=If[And@@(StringQ/@files), FFLoadEvaluationsImplem[GetGraphId[gid],files], $Failed];
FFSamplesFileSize[file_String]:=FFSamplesFileSizeImplem[file];
FFNParsFromDegreesFile[file_String]:=If[TrueQ[#==$Failed],$Failed,{"NParsIn"->#[[1]],"NParsOut"->#[[2]]}]&[FFNParsFromDegreesFileImpl[file]];


RegisterSimpleSubgraph[gid_,inputs_,{subgraphid_}]:=Catch[FFSimpleSubGraphImplem[gid,inputs,GetGraphId[subgraphid]]];
FFAlgSimpleSubgraph[gid_,id_,inputs_List,subgraphid_]:=FFRegisterAlgorithm[RegisterSimpleSubgraph,gid,id,inputs,{subgraphid}];
FFAlgSimpleSubgraph[gid_,id_,inputs_List]:=FFAlgSimpleSubgraph[gid,id,inputs,id];


RegisterMemoizedSubgraph[gid_,inputs_,{subgraphid_}]:=Catch[FFMemoizedSubGraphImplem[gid,inputs,GetGraphId[subgraphid]]];
FFAlgMemoizedSubgraph[gid_,id_,inputs_List,subgraphid_]:=FFRegisterAlgorithm[RegisterMemoizedSubgraph,gid,id,inputs,{subgraphid}];
FFAlgMemoizedSubgraph[gid_,id_,inputs_List]:=FFAlgMemoizedSubgraph[gid,id,inputs,id];


RegisterSubgraphMap[gid_,inputs_,{subgraphid_}]:=Catch[FFSubGraphMapImplem[gid,inputs,GetGraphId[subgraphid]]];
FFAlgSubgraphMap[gid_,id_,inputs_List,subgraphid_]:=FFRegisterAlgorithm[RegisterSubgraphMap,gid,id,inputs,{subgraphid}];
FFAlgSubgraphMap[gid_,id_,inputs_List]:=FFAlgSubgraphMap[gid,id,inputs,id];


RegisterAlgLaurent[gid_,inputs_,{subgraphid_,order_,maxdeg_}]:=Catch[FFAlgLaurentImplem[gid,inputs,GetGraphId[subgraphid],CheckedInt32List[order],toFFInternalUnsignedFlag["MaxDegree",maxdeg]]];
Options[FFAlgLaurent]:={"MaxDegree"->Automatic};
FFAlgLaurent[gid_,id_,inputs_List,subgraphid_,order_Integer,opt:OptionsPattern[]]:=FFAlgLaurent[gid,id,inputs,subgraphid,{order},opt];
FFAlgLaurent[gid_,id_,inputs_List,subgraphid_,order_List,OptionsPattern[]]:=FFRegisterAlgorithm[RegisterAlgLaurent,gid,id,inputs,{subgraphid,order,OptionValue["MaxDegree"]}];

PartitionsWithLen[l_, p_]:=PartitionsWithLenImpl[l,Accumulate[p]];
PartitionsWithLenImpl[l_, p_]:=Inner[l[[#1;;#2]]&,Join[{0},Most[p]]+1,p,List];
FFLaurentSol[solin_,expvar_,info_]:=Module[{sol,len},
  len = Inner[Max[#2-#1+1,0]&,info[[1]],info[[2]],List];
  sol = PartitionsWithLen[solin,len];
  SeriesData[expvar,0,sol[[#]],info[[1,#]],info[[2,#]]+1,1]&/@Range[Length[sol]]
];


RegisterDenseSolver[gid_,inputs_,{params_,eqsin_,vars_,neededvarsin_,applyfun_}]:=Module[
  {eqs, neededvars, lincoeffs, dummy},
  Catch[
    eqs = Select[eqsin, !TrueQ[#]&];
    neededvars = If[TrueQ[neededvarsin==Automatic], vars, neededvarsin];
    CheckVariables[vars];
    CheckVariables[neededvars];
    If[!SubsetQ[vars,neededvars], Message[FF::badneededvars]; Throw[$Failed];];
    If[!TrueQ[eqs[[0]]==List && And@@(((#[[0]] == Equal) && (Length[#] == 2))&/@eqs)],
        Message[FF::badsystem]; Throw[$Failed];
    ];
    eqs = (#[[1]]-#[[2]])&/@eqs;
    lincoeffs = LinearEqCoeffs[#,vars,applyfun]&/@eqs;
    If[!TrueQ[params == {}],
      lincoeffs = Map[(toFFInternalRatFun[#,params]&/@#)&, lincoeffs];
      FFRegisterDenseSolverImplem[gid,inputs,Length[params],Length[vars],lincoeffs,((Position[vars,#][[1,1]])&/@neededvars)-1],
      lincoeffs = Map[(ToString[#,InputForm]&/@#)&, lincoeffs];
      FFRegisterDenseSolverImplemN[gid,inputs,Length[vars],lincoeffs,((Position[vars,#][[1,1]])&/@neededvars)-1]
    ]
  ]
];

Options[FFAlgDenseSolver]:={"NeededVars"->Automatic, "ApplyFunction"->Identity};
FFAlgDenseSolver[gid_,id_,inputs_List,params_,eqs_,vars_,OptionsPattern[]]:=Module[{},
  FFRegisterAlgorithm[RegisterDenseSolver, gid, id, inputs, {params, eqs, vars, OptionValue["NeededVars"], OptionValue["ApplyFunction"]}]
];


RegisterNodeDenseSolver[gid_,inputs_,{neqs_,vars_,neededvarsin_}]:=Module[
  {neededvars, lincoeffs, dummy},
  Catch[
    neededvars = If[TrueQ[neededvarsin==Automatic], vars, neededvarsin];
    CheckVariables[vars];
    CheckVariables[neededvars];
    If[!SubsetQ[vars,neededvars], Message[FF::badneededvars]; Throw[$Failed];];
    FFRegisterNodeDenseSolverImplem[gid,inputs,CheckedInt[neqs],Length[vars],((Position[vars,#][[1,1]])&/@neededvars)-1]
  ]
];

Options[FFAlgNodeDenseSolver]:={"NeededVars"->Automatic};
FFAlgNodeDenseSolver[gid_,id_,inputs_List,neqs_,vars_,OptionsPattern[]]:=Module[{},
  FFRegisterAlgorithm[RegisterNodeDenseSolver, gid, id, inputs, {neqs, vars, OptionValue["NeededVars"]}]
];


(*RegisterSparseSolver[gid_,inputs_,{params_,eqsin_,vars_,neededvarsin_,applyfun_}]:=Module[
  {eqs, neededvars, lincoeffs, dummy, position, tointernal},
  Catch[
    position = Association[{}];
    Table[position[vars[[ii]]]=ii-1;,{ii,Length[vars]}];
    eqs = Select[eqsin, !TrueQ[#]&];
    neededvars = If[TrueQ[neededvarsin==Automatic], vars, neededvarsin];
    CheckVariables[vars];
    CheckVariables[neededvars];
    If[!SubsetQ[vars,neededvars], Message[FF::badneededvars]; Throw[$Failed];];
    If[!TrueQ[eqs[[0]]==List && And@@(((#[[0]] == Equal) && (Length[#] == 2))&/@eqs)],
        Message[FF::badsystem]; Throw[$Failed];
    ];
    eqs = (#[[1]]-#[[2]])&/@eqs;
    lincoeffs = SparseLinearEqCoeffs[#,vars,applyfun]&/@eqs;
    If[!TrueQ[params == {}],
      tointernal[expr_]:=tointernal[expr]=toFFInternalRatFun[expr,params];,
      tointernal[expr_]:=tointernal[expr]=If[!TrueQ[FFRationalQ[expr]],
                                             Message[FF::badfunarg, expr, params]; Throw[$Failed],
                                             ToString[expr,InputForm]];
    ];
    lincoeffs = Map[{#[[1]],tointernal/@#[[2]]}&, lincoeffs];
    If[!TrueQ[params == {}],
      FFRegisterSparseSolverImplem[gid,inputs,Length[params],Length[vars],lincoeffs,position/@neededvars],
      FFRegisterSparseSolverImplemN[gid,inputs,Length[vars],lincoeffs,position/@neededvars]
    ]
  ]
];*)

RegisterSparseSolver[gid_,inputs_,{params_,eqsin_,vars_,neededvarsin_,applyfun_}]:=Module[
  {varmap,xx,xvars},
  Catch[
    xvars = xx/@Range[Length[vars]];
    varmap = Dispatch[Inner[Rule, vars, xvars, List]];
    RegisterSparseSolver[gid,inputs, {params,eqsin,xvars,(Union[Cases[{#},_xx,Infinity]]&),neededvarsin,applyfun}/.varmap]
  ]
];

RegisterSparseSolver[gid_,inputs_,{params_,eqsin_,vars_,pattern_,neededvarsin_,applyfun_}]:=Module[
  {eqs, neededvars, lincoeffs, dummy, position, tointernal},
  Catch[
    position = Association[{}];
    Table[position[vars[[ii]]]=ii-1;,{ii,Length[vars]}];
    eqs = Select[eqsin, !TrueQ[#]&];
    neededvars = If[TrueQ[neededvarsin==Automatic], vars, neededvarsin];
    CheckVariables[vars];
    CheckVariables[neededvars];
    If[!SubsetQ[vars,neededvars], Message[FF::badneededvars]; Throw[$Failed];];
    If[!TrueQ[eqs[[0]]==List && And@@(((#[[0]] == Equal) && (Length[#] == 2))&/@eqs)],
        Message[FF::badsystem]; Throw[$Failed];
    ];
    eqs = (#[[1]]-#[[2]])&/@eqs;
    lincoeffs = SparseLinearEqCoeffsWithPattern[#,Length[vars],pattern,position,applyfun]&/@eqs;
    If[!FreeQ[lincoeffs,Missing["KeyAbsent",_]],
        Message[FF::badpattern,#[[2]]&/@Union[Cases[lincoeffs,_Missing,Infinity]]]; Throw[$Failed];
    ];
    If[!TrueQ[params == {}],
      tointernal[expr_]:=tointernal[expr]=toFFInternalRatFun[expr,params];,
      tointernal[expr_]:=tointernal[expr]=If[!TrueQ[FFRationalQ[expr]],
                                             Message[FF::badfunarg, expr, params]; Throw[$Failed],
                                             ToString[expr,InputForm]];
    ];
    lincoeffs = Map[{#[[1]],tointernal/@#[[2]]}&, lincoeffs];
    If[!TrueQ[params == {}],
      FFRegisterSparseSolverImplem[gid,inputs,Length[params],Length[vars],lincoeffs,position/@neededvars],
      FFRegisterSparseSolverImplemN[gid,inputs,Length[vars],lincoeffs,position/@neededvars]
    ]
  ]
];

Options[FFAlgSparseSolver]:={"NeededVars"->Automatic, "ApplyFunction"->Identity, "VarsPattern"->Automatic};
FFAlgSparseSolver[gid_,id_,inputs_List,params_,eqs_,vars_,OptionsPattern[]]:=Module[{},
  If[TrueQ[OptionValue["VarsPattern"]==Automatic],
   FFRegisterAlgorithm[RegisterSparseSolver, gid, id, inputs, {params, eqs, vars, OptionValue["NeededVars"], OptionValue["ApplyFunction"]}],
   FFRegisterAlgorithm[RegisterSparseSolver, gid, id, inputs, {params, eqs, vars, OptionValue["VarsPattern"], OptionValue["NeededVars"], OptionValue["ApplyFunction"]}]
  ]
];


Options[FFSerializeSparseEqs]={"ApplyFunction"->Identity};
FFSerializeSparseEqs[filename_,params_,eqsin_,vars_,pattern_,position_,OptionsPattern[]]:=Module[
  {eqs, lincoeffs, dummy, tointernal,applyfun},
  Catch[
    applyfun=OptionValue["ApplyFunction"];
    eqs = Select[eqsin, !TrueQ[#]&];
    CheckVariables[vars];
    If[!TrueQ[eqs[[0]]==List && And@@(((#[[0]] == Equal) && (Length[#] == 2))&/@eqs)],
        Message[FF::badsystem]; Throw[$Failed];
    ];
    eqs = (#[[1]]-#[[2]])&/@eqs;
    lincoeffs = SparseLinearEqCoeffsWithPattern[#,Length[vars],pattern,position,applyfun]&/@eqs;
    If[!FreeQ[lincoeffs,Missing["KeyAbsent",_]],
        Message[FF::badpattern,#[[2]]&/@Union[Cases[lincoeffs,_Missing,Infinity]]]; Throw[$Failed];
    ];
    tointernal[expr_]:=tointernal[expr]=toFFInternalRatFun[expr,params];
    lincoeffs = Map[{#[[1]],tointernal/@#[[2]]}&, lincoeffs];
    Export[filename,lincoeffs,"MX"]
  ]
];

RegisterSparseSerializedEqs[gid_, inputs_, {params_, files_, vars_, neededvarsin_}]:=Module[
  {needed,position,neededvars,lincoeffs},
  position = Association[{}];
  Table[position[vars[[ii]]]=ii-1;,{ii,Length[vars]}];
  neededvars = If[TrueQ[neededvarsin==Automatic], vars, neededvarsin];
  CheckVariables[vars];
  CheckVariables[neededvars];
  If[!SubsetQ[vars,neededvars], Message[FF::badneededvars]; Throw[$Failed];];
  lincoeffs=Join@@(Import/@files);
  FFRegisterSparseSolverImplem[gid,inputs,Length[params],Length[vars],lincoeffs,position/@neededvars]
];

Options[FFAlgSerializedSparseSolver]:={"NeededVars"->Automatic};
FFAlgSerializedSparseSolver[gid_,id_,inputs_List,params_,files_List,vars_,OptionsPattern[]]:=FFRegisterAlgorithm[RegisterSparseSerializedEqs, gid, id, inputs, {params, files, vars, OptionValue["NeededVars"]}];


Options[FFSparseEqsToJSON]={"ApplyFunction"->Identity};
FFSparseEqsToJSON[filename_,params_,eqsin_,vars_,pattern_,position_,OptionsPattern[]]:=Module[
  {eqs, lincoeffs, dummy, tointernal,applyfun,nvars},
  Catch[
    applyfun=OptionValue["ApplyFunction"];
    eqs = Select[eqsin, !TrueQ[#]&];
    If[vars[[0]] == List,
      CheckVariables[vars];
      nvars = Length[vars];,
      nvars = vars;
    ];
    If[!TrueQ[eqs[[0]]==List && And@@(((#[[0]] == Equal) && (Length[#] == 2))&/@eqs)],
        Message[FF::badsystem]; Throw[$Failed];
    ];
    eqs = (#[[1]]-#[[2]])&/@eqs;
    lincoeffs = SparseLinearEqCoeffsWithPattern[#,nvars,pattern,position,applyfun]&/@eqs;
    If[!FreeQ[lincoeffs,Missing["KeyAbsent",_]],
        Message[FF::badpattern,#[[2]]&/@Union[Cases[lincoeffs,_Missing,Infinity]]]; Throw[$Failed];
    ];
    tointernal[expr_]:=tointernal[expr]=toFFJSONRatFun[expr,params];
    lincoeffs = Map[{Length[#[[1]]],#[[1]],tointernal/@#[[2]]}&, lincoeffs];
    Export[filename,{Length[lincoeffs],lincoeffs},"RawJSON","Compact"->True]
  ]
];

Options[FFSparseSystemToJSON]:={"NeededVars"->Automatic};
FFSparseSystemToJSON[file_,neqs_,vars_,pars_,files_List,OptionsPattern[]]:=Module[
  {neededlist,position},
  position = Association[{}];
  Table[position[vars[[ii]]]=ii-1;,{ii,Length[vars]}];
  neededlist = OptionValue["NeededVars"];
  If[neededlist==Automatic, neededlist = vars];
  CheckVariables[vars];
  CheckVariables[neededlist];
  If[!SubsetQ[vars,neededlist], Message[FF::badneededvars]; Throw[$Failed];];
  Export[file,{neqs,Length[vars],Length[pars],Length[neededlist],position/@neededlist,Length[files],files},"RawJSON","Compact"->True]
];

RegisterSparseEqsFromJSON[gid_,inputs_,{file_}]:=FFAlgJSONSparseSolverImplem[gid,inputs,file];

Options[FFAlgJSONSparseSolver]:={"NeededVars"->Automatic};
FFAlgJSONSparseSolver[gid_,id_,inputs_List,file_,OptionsPattern[]]:=FFRegisterAlgorithm[RegisterSparseEqsFromJSON, gid, id, inputs, {file}];


FFRatFunToJSON[filename_,vars_,funcs_]:=Catch[Export[filename,{Length[funcs],Length[vars],toFFJSONRatFun[#,vars]&/@funcs},"RawJSON","Compact"->True]];
RegisterJSONRatFunEval[gid_,inputs_,{file_}]:=FFAlgJSONRatFunEvalImplem[gid,inputs,file];
FFAlgJSONRatFunEval[gid_,id_,inputs_List,file_]:=FFRegisterAlgorithm[RegisterJSONRatFunEval, gid, id, inputs, {file}];


FFPolyFromJSON[expr_,vars_]:=Plus@@((ToExpression[#[[1]]] Times@@(vars^#[[2]]))&/@expr);
FFRatFunFromJSONImpl[expr_,vars_]:=FFPolyFromJSON[expr[[1,2]],vars]/FFPolyFromJSON[expr[[2,2]],vars];
FFRatFunFromJSONImpl[0,vars_]:=0;
FFRatFunFromJSON[file_,vars_]:=FFRatFunFromJSONImpl[#,vars]&/@Import[file][[3]];


FFUIntListToJSON[file_,list_]:=Catch[Export[file,{Length[list],CheckedUInt32List[list]},"RawJSON","Compact"->True]];


FFSolverResetNeededVars[gid_,id_,vars_,needed_]:=Module[{position},
  position = Association[{}];
  Table[position[vars[[ii]]]=ii-1;,{ii,Length[vars]}];
  FFSolverResetNeededVarsImplem[GetGraphId[gid],GetAlgId[gid,id],position/@needed]
];


FFListSubsetToJSON[file_,list_,subset_]:=Module[{position,ii},
  position=Association[{}];
  Do[position[list[[ii]]]=ii-1;,{ii,Length[list]}];
  FFUIntListToJSON[file,position/@subset]
];
FFListSubsetFromJSON[file_,list_]:=list[[Import[file][[2]]+1]];


FFSolverOnlyHomogeneous[gid_,id_]:=Catch[FFSolverOnlyHomogeneousImplem[GetGraphId[gid],GetAlgId[gid,id]]];


FFSolverSparseOutput[gid_,id_]:=Catch[FFSolverSparseOutputImplem[GetGraphId[gid],GetAlgId[gid,id]]];


FFLearn[gid_]:=Catch[FFLearnImplem[GetGraphId[gid]]];


FFLaurentLearn[graph_]:=FFLearn[graph];


FFDenseSolverLearn[gid_,vars_]:=Module[
  {depv,indepv,zerov,learn},
  Catch[
    learn = FFLearnImplem[GetGraphId[gid]];
    If[!TrueQ[learn[[0]]==List], Throw[learn]];
    {depv,indepv,zerov} = learn;
    If[Length[depv]==0 && Length[indepv]==0  && Length[zerov]==0, Return[FFImpossible]];
    {"DepVars"->vars[[depv+1]],"IndepVars"->vars[[indepv+1]],"ZeroVars"->vars[[zerov+1]]}
  ]
];


FFSparseSolverLearn[gid_,vars_]:=Module[
  {depv,indepv,learn,sparseout,varswc},
  Catch[
    learn = FFLearnImplem[GetGraphId[gid]];
    If[!TrueQ[learn[[0]]==List], Throw[learn]];
    {depv,indepv,sparseout} = learn;
    If[Length[depv]==0 && Length[indepv]==0, Return[FFImpossible]];
    If[!TrueQ[sparseout==1],
      {"DepVars"->vars[[depv+1]],"IndepVars"->vars[[indepv+1]],"SparseOutput"->False},
      varswc=Join[vars,{1}];
      {"DepVars"->vars[[depv+1]],"IndepVars"->(varswc[[#]]&/@(indepv+1)),"SparseOutput"->True}
    ]
  ]
];


ConvertDenseLearn[learn_,vars_]:=Module[{depv,indepv,zerov},
    {depv,indepv,zerov} = learn;
    If[Length[depv]==0 && Length[indepv]==0  && Length[zerov]==0, Return[FFImpossible]];
    {"DepVars"->vars[[depv+1]],"IndepVars"->vars[[indepv+1]],"ZeroVars"->vars[[zerov+1]]}
];


FFMultiFitLearn[gid_,vars_]:=Module[
  {learn},
  learn = FFLearn[gid];
  If[!TrueQ[learn[[0]]==List], Return[learn]];
  MapThread[ConvertDenseLearn,{learn,vars},1]
];


FFDenseSolverGetInfo[gid_,id_,vars_]:=Module[
  {depv,indepv,zerov},
  Catch[
    {depv,indepv,zerov} = FFAlgorithmGetInfoImplem[GetGraphId[gid],GetAlgId[gid,id]];
    If[Length[depv]==0 && Length[indepv]==0  && Length[zerov]==0, Return[FFImpossible]];
    {"DepVars"->vars[[depv+1]],"IndepVars"->vars[[indepv+1]],"ZeroVars"->vars[[zerov+1]]}
  ]
];


FFSparseSolverGetInfo[gid_,id_,vars_]:=Module[
  {depv,indepv,sparseout,varswc},
  Catch[
    {depv,indepv,sparseout} = FFAlgorithmGetInfoImplem[GetGraphId[gid],GetAlgId[gid,id]];
    If[Length[depv]==0 && Length[indepv]==0, Return[FFImpossible]];
    If[!TrueQ[sparseout==1],
      {"DepVars"->vars[[depv+1]],"IndepVars"->vars[[indepv+1]],"SparseOutput"->False},
      varswc=Join[vars,{1}];
      {"DepVars"->vars[[depv+1]],"IndepVars"->(varswc[[#]]&/@(indepv+1)),"SparseOutput"->True}
    ]
  ]
];


FFMultiFitGetInfo[gid_,id_,vars_]:=Module[
  {info},
  info = FFAlgorithmGetInfoImplem[GetGraphId[gid],GetAlgId[gid,id]];
  MapThread[ConvertDenseLearn,{info,vars},1]
];


FFSubgraphReconstructGetInfo[gid_,id_,vars_]:=Module[
  {info},
  info = FFAlgorithmGetInfoImplem[GetGraphId[gid],GetAlgId[gid,id]];
  (((Times@@(vars^#))&/@#)&/@#)&/@info
];


FFSolverNIndepEqs[gid_,id_]:=FFSolverNIndepEqsImplem[GetGraphId[gid],GetAlgId[gid,id]];


FFSolverIndepEqs[gid_,id_]:=If[TrueQ[#[[0]]==List],#+1,#]&[FFSolverIndepEqsImplem[GetGraphId[gid],GetAlgId[gid,id]]];


FFSparseSolverMarkAndSweepEqs[gid_,id_]:=FFSparseSolverMarkAndSweepEqsImplem[GetGraphId[gid],GetAlgId[gid,id]];


FFSparseSolverDeleteUnneededEqs[gid_,id_]:=FFSparseSolverDeleteUnneededEqsImplem[GetGraphId[gid],GetAlgId[gid,id]];


RegisterAlgRatFunEval[gid_,inputs_,{vars_List,functions_}]:=Catch[FFAlgRatFunEvalImplem[gid,inputs,Length[vars],toFFInternalRatFun[#,vars]&/@functions]];
FFAlgRatFunEval[gid_,id_,inputs_List,params_,functions_List]:=FFRegisterAlgorithm[RegisterAlgRatFunEval,gid,id,inputs,{params,functions}];


RegisterAlgRatFunEvalFromCoeffs[gid_,inputs_,{coeffs_,vars_List,functions_}]:=Module[
  {map},
  map = Association[{}];
  Do[map[coeffs[[ii]]] = ii-1;,{ii,Length[coeffs]}];
  Catch[
    FFAlgRatFunEvalFromCoeffsImplem[gid,inputs,Length[vars],Length[coeffs],toFFInternalRatFunCoeffMap[#,vars,map]&/@functions]
  ]
];
FFAlgRatFunEvalFromCoeffs[gid_,id_,inputs_List,coeffs_,params_,functions_List]:=FFRegisterAlgorithm[RegisterAlgRatFunEvalFromCoeffs,gid,id,inputs,{coeffs,params,functions}];


RegisterAlgRatNumEval[gid_,{},{numbers_}]:=Catch[If[!AllTrue[numbers,FFRationalQ],Throw[$Failed]]; FFAlgRatNumEvalImplem[gid,{},ToString[#,InputForm]&/@numbers]];
FFAlgRatNumEval[gid_,id_,numbers_List]:=FFRegisterAlgorithm[RegisterAlgRatNumEval,gid,id,{},{numbers}];


RegisterAlgChain[gid_,inputs_,{}]:=Catch[FFAlgChainImplem[gid,inputs]];
FFAlgChain[gid_,id_,inputs_List]:=FFRegisterAlgorithm[RegisterAlgChain,gid,id,inputs,{}];


ValidateTakeElemsList[a_List]:=If[AllTrue[a,#[[0]]==List && Length[#]==2&],a,Throw[$Failed]];
TakeElemsToInternal[a_List]:=ValidateTakeElemsList[a]-1;
TakeElemsToInternal[full_List->subset_List]:=Module[{position,i,j,sublist,res},
  If[!AllTrue[full,#[[0]]==List&],Throw[$Failed]];
  position=Association[{}];
  Do[
    sublist = full[[i]];
    Do[position[sublist[[j]]] = {i,j};,{j,Length[sublist]}];
  ,{i,Length[full]}];
  res=position/@subset;
  If[!FreeQ[res,Missing],Throw[$Failed]];
  res-1
];
RegisterAlgTake[gid_,inputs_,{elems_}]:=Catch[FFAlgTakeImplem[gid,inputs,Flatten[TakeElemsToInternal[elems]]]];
FFAlgTake[gid_,id_,inputs_List,elems_]:=FFRegisterAlgorithm[RegisterAlgTake,gid,id,inputs,{elems}];


MultiTakeElemsToInternal[a_List]:=Flatten[TakeElemsToInternal[#]]&/@a;
MultiTakeElemsToInternal[full_List->subsets_List]:=Flatten[TakeElemsToInternal[full->#]]&/@subsets;
RegisterAlgTakeAndAdd[gid_,inputs_,{elems_}]:=Catch[FFAlgTakeAndAddImplem[gid,inputs,MultiTakeElemsToInternal[elems]]];
FFAlgTakeAndAdd[gid_,id_,inputs_List,elems_]:=FFRegisterAlgorithm[RegisterAlgTakeAndAdd,gid,id,inputs,{elems}];


RegisterAlgSlice[gid_,inputs_,{start_,end_}]:=Catch[FFAlgSliceImplem[gid,inputs,start-1,end]];
FFAlgSlice[gid_,id_,inputs_List,start_,end_:-1]:=FFRegisterAlgorithm[RegisterAlgSlice,gid,id,inputs,{start,end}];


RegisterAlgAdd[gid_,inputs_,{}]:=Catch[FFAlgAddImplem[gid,inputs]];
FFAlgAdd[gid_,id_,inputs_List]:=FFRegisterAlgorithm[RegisterAlgAdd,gid,id,inputs,{}];


RegisterAlgMul[gid_,inputs_,{}]:=Catch[FFAlgMulImplem[gid,inputs]];
FFAlgMul[gid_,id_,inputs_List]:=FFRegisterAlgorithm[RegisterAlgMul,gid,id,inputs,{}];


RegisterAlgMatMul[gid_,inputs_,{rows1_,cols1_,cols2_}]:=Catch[FFAlgMatMulImplem[gid,inputs,CheckedInt32[rows1],CheckedInt32[cols1],CheckedInt32[cols2]]];
FFAlgMatMul[gid_,id_,inputs_List,rows1_,cols1_,cols2_]:=FFRegisterAlgorithm[RegisterAlgMatMul,gid,id,inputs,{rows1,cols1,cols2}];


ValidateSparseMatElemsList[a_List]:=If[AllTrue[a,#[[0]]==List&],a,Throw[$Failed]];
GetSparseColumns[a_List] := ValidateSparseMatElemsList[a]-1;
GetSparseColumns[full_->columns_] := #[[2]]&/@TakeElemsToInternal[{full}->#]&/@columns;
RegisterAlgSparseMatMul[gid_,inputs_,{rows1_,cols1_,cols2_,columns1_,columns2_}]:=Catch[FFAlgSparseMatMulImplem[gid,inputs,CheckedInt32[rows1],CheckedInt32[cols1],CheckedInt32[cols2],GetSparseColumns[columns1],GetSparseColumns[columns2]]];
FFAlgSparseMatMul[gid_,id_,inputs_List,rows1_,cols1_,cols2_,columns1_,columns2_]:=FFRegisterAlgorithm[RegisterAlgSparseMatMul,gid,id,inputs,{rows1,cols1,cols2,columns1,columns2}];


RegisterAlgNonZeroes[gid_,inputs_,{}]:=FFAlgNonZeroesImplem[gid,inputs];
FFAlgNonZeroes[gid_,id_,inputs_List]:=FFRegisterAlgorithm[RegisterAlgNonZeroes,gid,id,inputs,{}];


FFNonZeroesGetInfo[gid_,id_]:={"All"->#[[1]],"NonZero"->(#[[2]]+1)}&[FFAlgorithmGetInfoImplem[GetGraphId[gid],GetAlgId[gid,id]]];
FFNonZeroesLearn[gid_]:=If[Length[#]==0,#,{"All"->#[[1]],"NonZero"->(#[[2]]+1)}]&[FFLearnImplem[GetGraphId[gid]]];


FFIndependentOf[id_, vars_List, var_]:=FFIndependentOfImplem[GetGraphId[id],Position[vars,var][[1,1]]-1];


Options[FFReconstructFunction]:=Join[Options[FFAlgorithmSetReconstructionOptions],{"NThreads"->FFNThreads,"MinPrimes"->1}];
FFReconstructFunction[id_,vars_,OptionsPattern[]] := Module[
  {np,maxnp,opt,res,nthreads,tmp,thisopt},
  opt = (#[[1]]->OptionValue[#[[1]]])&/@Options[FFReconstructFunction];
  maxnp = OptionValue["MaxPrimes"];
  If[TrueQ[maxnp==Automatic], maxnp = 5];
  nthreads = OptionValue["NThreads"];
  If[TrueQ[nthreads==Automatic], nthreads = If[TrueQ[Length[vars]==1],1,FFAutomaticNThreads[]]];
  If[Length[vars]==1,
    thisopt = Join[{"MaxPrimes"->maxnp},FilterRules[opt,Select[Options[FFReconstructUnivariate],FreeQ[#,"MaxPrimes"]&]]];
    Return[FFReconstructUnivariate[id,vars,Sequence@@thisopt]]
  ];
  res = FFAllDegrees[id, nthreads, Sequence@@FilterRules[opt,Options[FFAllDegrees]]];
  If[TrueQ[res == $Failed], Return[$Failed]];
  np = OptionValue["MinPrimes"];
  tmp = FFMissingPoints;
  res = Catch[While[True,
    If[TrueQ[np>maxnp], Throw[tmp]];
    thisopt = Join[{"MaxPrimes"->np},FilterRules[opt,Select[Options[FFSample],FreeQ[#,"MaxPrimes"]&]]];
    FFSample[id,nthreads,Sequence@@thisopt];
    tmp = FFReconstructFromCurrentEvaluations[id,vars,nthreads,Sequence@@thisopt];
    If[!TrueQ[tmp[[0]]==List],
      Switch[tmp,
             FFMissingPrimes, np=np+1;,
             _, Throw[tmp];
      ];,
      Throw[tmp];
   ]
  ]];
  res
];


Options[FFReconstructFunctionMod]:=Join[Options[FFAlgorithmSetReconstructionOptions],{"NThreads"->FFNThreads}];
FFReconstructFunctionMod[id_,vars_,OptionsPattern[]] := Module[
  {np,maxnp,opt,res,nthreads,tmp,thisopt},
  opt = (#[[1]]->OptionValue[#[[1]]])&/@Options[FFReconstructFunctionMod];
  maxnp = 1;
  If[TrueQ[maxnp==Automatic], maxnp = 5];
  nthreads = OptionValue["NThreads"];
  If[TrueQ[nthreads==Automatic], nthreads = If[TrueQ[Length[vars]==1],1,FFAutomaticNThreads[]]];
  If[Length[vars]==1,
    thisopt = Join[{"MaxPrimes"->maxnp},FilterRules[opt,Select[Options[FFReconstructUnivariate],FreeQ[#,"MaxPrimes"]&]]];
    Return[FFReconstructUnivariateMod[id,vars,Sequence@@thisopt]]
  ];
  res = FFAllDegrees[id, nthreads, Sequence@@FilterRules[opt,Options[FFAllDegrees]]];
  If[TrueQ[res == $Failed], Return[$Failed]];
  np = 1;
  tmp = FFMissingPoints;
  res = Catch[
    thisopt = Join[{"MaxPrimes"->np},FilterRules[opt,Select[Options[FFSample],FreeQ[#,"MaxPrimes"]&]]];
    FFSample[id,nthreads,Sequence@@thisopt];
    tmp = FFReconstructFromCurrentEvaluationsMod[id,vars,nthreads,Sequence@@thisopt];
    Throw[tmp];
  ];
  res
];


Options[FFParallelReconstructUnivariate]=Join[Options[FFReconstructFunction],{"MinDegree"->Automatic,"DegreeStep"->Automatic}];
FFParallelReconstructUnivariate[id_,vars_,OptionsPattern[]]:=Module[
  {opt,mindeg,maxdeg,degstep,deg,maxsp,maxnp,nthreads,sp,np,tmp,res,thisopt},
  opt = (#[[1]]->OptionValue[#[[1]]])&/@Options[FFParallelReconstructUnivariate];
  maxsp = OptionValue["MaxSingularPoints"];
  If[TrueQ[maxsp==Automatic], maxsp = 1];
  maxnp = OptionValue["MaxPrimes"];
  If[TrueQ[maxnp==Automatic], maxnp = 5];
  nthreads = OptionValue["NThreads"];
  If[TrueQ[nthreads==Automatic], nthreads = FFAutomaticNThreads[]];
  mindeg = OptionValue["MinDegree"];
  If[TrueQ[mindeg == Automatic], mindeg = nthreads];
  maxdeg = OptionValue["MaxDegree"];
  If[TrueQ[maxdeg == Automatic], maxdeg = FFDefaultMaxDegImplem[]];
  degstep = OptionValue["DegreeStep"];
  If[TrueQ[degstep == Automatic], degstep = nthreads];
  sp = 0;
  np = OptionValue["MinPrimes"];
  deg = mindeg;
  tmp = FFMissingPoints;
  res = Catch[While[True,
    If[TrueQ[sp>maxsp || np>maxnp || deg > maxdeg], Throw[tmp]];
    thisopt = Join[{"MaxSingularPoints"->sp,"MaxPrimes"->np,"MaxDegree"->deg},FilterRules[opt,Select[Options[FFSample],FreeQ[#,"MaxSingularPoints"|"MaxPrimes"|"MaxDegree"]&]]];
    FFSample[id,nthreads,Sequence@@thisopt];
    tmp = FFReconstructFromCurrentEvaluations[id,vars,nthreads,Sequence@@thisopt];
    If[!TrueQ[tmp[[0]]==List],
      Switch[tmp,
             FFMissingPrimes, np=np+1;,
             FFMissingPoints, sp=sp+1;,
             $Failed, If[deg<maxdeg, deg = deg+degstep; If[deg>maxdeg,deg=maxdeg];, deg = deg+degstep;];,
             _, Throw[tmp];
      ];,
      Throw[tmp];
   ]
  ]];
  res
];


SparseSolWSparseOut[sol_,depv_,indepv_]:=Module[{psol,rhs},
  psol=PartitionsWithLen[sol,Length/@indepv];
  rhs=Table[psol[[ii]].indepv[[ii]],{ii,Length[psol]}];
  Inner[Rule,depv,rhs,List]
];


FFDenseSolverSol[sol_,learninfo_]:=Module[{depv,indepv,zero},
  {depv,indepv,zero}={"DepVars","IndepVars","ZeroVars"}/.learninfo;
  If[!TrueQ[Length[depv]==0],
    Join[Inner[Rule,depv,((#).Join[indepv,{1}])&/@ArrayReshape[sol,{Length[depv],Length[indepv]+1}],List],(#->0)&/@zero],
    (#->0)&/@zero
  ]
];
FFSparseSolverSol[sol_,learninfo_]:=Module[{depv,indepv,sparseout},
  {depv,indepv,sparseout}={"DepVars","IndepVars","SparseOutput"}/.learninfo;
  If[TrueQ[sparseout],Return[SparseSolWSparseOut[sol,depv,indepv]]];
  If[!TrueQ[Length[depv]==0],
    Inner[Rule,depv,((#).Join[indepv,{1}])&/@ArrayReshape[sol,{Length[depv],Length[indepv]+1}],List],
    {}
  ]
];
FFNonZeroesSol[sol_,learninfo_]:=Module[{tot,nonzero,ret},
  {tot,nonzero}={"All","NonZero"}/.learninfo;
  ret = SparseArray[{},{tot},0];
  ret[[nonzero]] = sol;
  ret
];
FFMultiFitSol[sol_,learn_]:=MapThread[FFDenseSolverSol,{PartitionsWithLen[sol,(Length["DepVars"/.#]*((Length["IndepVars"]/.#)+1))&/@learn],learn},1];


AutoReconstructionOptions[]:=Options[FFReconstructFunction];


Options[FFDenseSolve] := Join[{"Parameters"->Automatic, "IndepVarsOnly"->False},
                                   AutoReconstructionOptions[],
                                   Options[FFAlgDenseSolver]];
FFDenseSolve[eqs_, vars_, OptionsPattern[]] := Module[
	{params,res,graph,in,sys,learn,opt},
	opt = (#[[1]]->OptionValue[#[[1]]])&/@Options[FFDenseSolve];
	FFNewGraph[graph];
    res = Catch[
      params = OptionValue["Parameters"];
      params = If[TrueQ[params==Automatic],
                  Complement[Variables[({#[[1]],#[[2]]})&/@eqs], vars],
                  params];
      If[!TrueQ[params == {}], CheckVariables[params]];
      FFGraphInputVars[graph,in,params];
      
      res = FFAlgDenseSolver[graph,sys,{in},params,eqs,vars,
                                Sequence@@FilterRules[{opt}, Options[FFAlgDenseSolver]]];
      If[res==$Failed,Throw[$Failed]];
      
      FFGraphOutput[graph,sys];
      learn = FFDenseSolverLearn[graph,vars];
      If[!TrueQ[learn[[0]]==List],Throw[learn]];
      If[TrueQ[OptionValue["IndepVarsOnly"]], Throw["IndepVars"/.learn]];
      
      res = If[TrueQ[params == {}],
                FFReconstructNumeric[graph, Sequence@@FilterRules[{opt}, Options[FFReconstructNumeric]]],
                FFReconstructFunction[graph,params, Sequence@@FilterRules[{opt}, Options[FFReconstructFunction]]]
             ];
       If[!TrueQ[res[[0]]==List],Throw[res]];
       FFDenseSolverSol[res,learn]
    ];
    FFDeleteGraph[graph];
    res
];


Options[FFSparseSolve] := Join[{"Parameters"->Automatic, "IndepVarsOnly"->False, "MarkAndSweep"->True, "SparseOutput"->False},
                                   AutoReconstructionOptions[],
                                   Options[FFAlgSparseSolver]];
FFSparseSolve[eqs_, vars_, OptionsPattern[]] := Module[
	{params,res,graph,in,sys,learn,opt},
	opt = (#[[1]]->OptionValue[#[[1]]])&/@Options[FFSparseSolve];
	FFNewGraph[graph];
    res = Catch[
      params = OptionValue["Parameters"];
      params = If[TrueQ[params==Automatic],
                  Complement[Variables[({#[[1]],#[[2]]})&/@eqs], vars],
                  params];
      If[!TrueQ[params == {}], CheckVariables[params]];
      FFGraphInputVars[graph,in,params];
      
      res = FFAlgSparseSolver[graph,sys,{in},params,eqs,vars,
                                Sequence@@FilterRules[{opt}, Options[FFAlgSparseSolver]]];
      If[res==$Failed,Throw[$Failed]];
      If[(!TrueQ[OptionValue["IndepVarsOnly"]]) && TrueQ[OptionValue["SparseOutput"]],
        FFSolverSparseOutput[graph,sys];
      ];
      
      FFGraphOutput[graph,sys];
      learn = FFSparseSolverLearn[graph,vars];
      If[!TrueQ[learn[[0]]==List],Throw[learn]];
      If[TrueQ[OptionValue["IndepVarsOnly"]], Throw["IndepVars"/.learn]];
      If[TrueQ[OptionValue["MarkAndSweep"]],
        FFSparseSolverMarkAndSweepEqs[graph,sys];
        FFSparseSolverDeleteUnneededEqs[graph,sys];
      ];
      
      res = If[TrueQ[params == {}],
               FFReconstructNumeric[graph, Sequence@@FilterRules[{opt}, Options[FFReconstructNumeric]]],
               FFReconstructFunction[graph,params, Sequence@@FilterRules[{opt}, Options[FFReconstructFunction]]]
             ];
       If[!TrueQ[res[[0]]==List],Throw[res]];
       FFSparseSolverSol[res,learn]
    ];
    FFDeleteGraph[graph];
    res
];


Options[FFInverse] := FilterRules[Join[Options[FFDenseSolve],{"Sparse"->False}], Except["NeededVars"|"IndepVarsOnly"]];
FFInverse[mat_List, OptionsPattern[]]:=Module[
    {eqs, len, varx, vary, varsx, varsy, sol, params,res,graph,in,sys,learn,sparse,opt},
    opt = (#[[1]]->OptionValue[#[[1]]])&/@Options[FFInverse];
    
    len = Length[mat];
    If[!TrueQ[And@@((Length[#]==len)&/@mat)], Message[FF::badsquaremat]; Return[$Failed]];
    
    varsx = varx/@Range[len];
    varsy = vary/@Range[len];
    
    eqs = Table[mat[[ii]].varsx-vary[ii]==0,{ii,len}];
    
    sparse=OptionValue["Sparse"];
    
	FFNewGraph[graph];
    res = Catch[
      params = OptionValue["Parameters"];
      params = If[TrueQ[params==Automatic],
                  Variables[mat],
                  params];
      If[!TrueQ[params == {}], CheckVariables[params]];
      FFGraphInputVars[graph,in,params];
     
      res = If[sparse,
               FFAlgSparseSolver[graph,sys,{in},params,eqs,Join[varsx,varsy],
                                    "VarsPattern"->(DeleteDuplicates[Cases[{#},(_varx|_vary),Infinity]]&),
                                     Sequence@@FilterRules[{opt}, Options[FFAlgSparseSolver]]],
               FFAlgDenseSolver[graph,sys,{in},params,eqs,Join[varsx,varsy],
                                   Sequence@@FilterRules[{opt}, Options[FFAlgDenseSolver]]]];
      If[res==$Failed,Throw[$Failed]];
      FFSolverOnlyHomogeneous[graph,sys];
      
      FFGraphOutput[graph,sys];
      learn = If[sparse,
                 FFSparseSolverLearn[graph,Join[varsx,varsy]],
                 FFDenseSolverLearn[graph,Join[varsx,varsy]]];
      If[TrueQ[learn==FFImpossible],Throw[FFSingularMatrix]];
      If[!TrueQ[learn[[0]]==List],Throw[learn]];
      If[!TrueQ[("DepVars"/.learn) == varsx && ("IndepVars"/.learn) == varsy],Throw[FFSingularMatrix]];
      
      res = If[TrueQ[params == {}],
               FFReconstructNumeric[graph, Sequence@@FilterRules[{opt}, Options[FFReconstructNumeric]]],
               FFReconstructFunction[graph,params, Sequence@@FilterRules[{opt}, Options[FFReconstructFunction]]]
             ];
       If[!TrueQ[res[[0]]==List],Throw[res]];
       ArrayReshape[res,{len,len}]
    ];
    FFDeleteGraph[graph];
    res
];


Options[FFGraphEvaluate]={"PrimeNo"->0};
FFGraphEvaluate[g_,x_,OptionsPattern[]]:=FFGraphEvaluateImplem[GetGraphId[g],CheckedInt64List[x],CheckedInt32[OptionValue["PrimeNo"]]];


Options[FFGraphEvaluateMany]={"PrimeNo"->0,"NThreads"->FFNThreads};
FFGraphEvaluateMany[g_,x_List,OptionsPattern[]]:=FFGraphEvaluateListImplem[GetGraphId[g],If[TrueQ[#==Automatic],FFAutomaticNThreads[],CheckedInt32[#]]&@OptionValue["NThreads"],CheckedInt32[OptionValue["PrimeNo"]],CheckedInt64List/@x];


FFPrimeNo[i_]:=FFPrimeNoImplem[CheckedInt32[i]];


NoEmptyList[a_]:=a;
NoEmptyList[{}]={1};


RegisterAlgLinearFit[gid_,inputs_,{params_,
                                 delta_,integrandin_,
                                 tauvarsin_,vars_,
                                 extraeqs_,lsubst_,
                                 neededvarsin_,applyfun_}] := Module[
    {tauvars, lincoeffs, lininternal, internalres, boolflags, uintflags, paramsin,
     dummysym, loopvars, lvinternal, tausubset, nloops,
     ni, coeffs, integrand, integrinternal, position, neededvars,
     integrlv, deltalv, systype, depvars, indepvars, ignoreddep, indepvarvec},
    Catch[
        (* this is roughly the same as the first part of CutSolve *)
        position = Association[{}];
        Table[position[vars[[ii]]]=ii-1;,{ii,Length[vars]}];
        neededvars = If[TrueQ[neededvarsin==Automatic], vars, neededvarsin];
        CheckVariables[vars];
        CheckVariables[neededvars];
        If[!SubsetQ[vars,neededvars], Message[FF::badneededvars]; Throw[$Failed];];
        tauvars = If[tauvarsin == {}, {dummysym}, tauvarsin];
        CheckVariables[tauvars];
        loopvars = #[[1]]&/@lsubst;
        If[!TrueQ[Intersection[loopvars,tauvars]=={}], Message[FF::badvarule]; Throw[$Failed]];
        integrand = Together/@If[TrueQ[integrandin[[0]]==List], Flatten[integrandin], {integrandin}];
        tausubset = Select[tauvars,!FreeQ[{delta,integrand},#]&];
        loopvars = If[TrueQ[loopvars == {}], tauvars, Join[loopvars,tausubset]];
        CheckVariables[loopvars];

        (* external spinor stuff *)
        If[!TrueQ[params=={}],CheckVariables[params]];

        (* system data *)
        Which[TrueQ[delta[[0]]==List && Length[delta]==Length[vars]+1],
          lincoeffs = applyfun/@delta;,
          TrueQ[delta[[0]]==List && Length[delta]==Length[vars]],
          lincoeffs = applyfun/@Join[delta,{0}];,
          True,
          lincoeffs = Together/@LinearEqCoeffs[#,vars,applyfun]&[delta];
        ];
        If[!TrueQ[lsubst=={}],
          If[!TrueQ[params=={}],
            lvinternal = Map[toFFInternalPolyRatFun[#,tauvars,params]&, Join[(#[[2]])&/@lsubst,tausubset]];,
            lvinternal = Map[toFFInternalRatFun[#,tauvars]&, Join[(#[[2]])&/@lsubst,tausubset]];
          ],
          lvinternal={}
        ];
        deltalv = Table[CheckedInt64List[NoEmptyList[Select[Range[Length[loopvars]], !FreeQ[term,loopvars[[#]]]&]]],{term,lincoeffs}];
        integrlv = Table[CheckedInt64List[NoEmptyList[Select[Range[Length[loopvars]], !FreeQ[integr,loopvars[[#]]]&]]],{integr,integrand}];
		If[!TrueQ[params=={}],
		  
		  lininternal = Map[toFFInternalPolyRatFun[#[[1]],loopvars[[#[[2]]]],params]&, Table[{lincoeffs[[i]],deltalv[[i]]},{i,Length[lincoeffs]}]];
          integrinternal = Map[toFFInternalPolyRatFun[#[[1]],loopvars[[#[[2]]]],params]&, Table[{integrand[[i]],integrlv[[i]]},{i,Length[integrand]}]];,

          lininternal = Map[toFFInternalRatFun[#[[1]],loopvars[[#[[2]]]]]&, Table[{lincoeffs[[i]],deltalv[[i]]},{i,Length[lincoeffs]}]];
          integrinternal = Map[toFFInternalRatFun[#[[1]],loopvars[[#[[2]]]]]&, Table[{integrand[[i]],integrlv[[i]]},{i,Length[integrand]}]];
        ];

        (* finally *)
        If[!TrueQ[params=={}],
          FFAlgLinearFitImplem[gid,inputs,Length[tauvars],lininternal,integrinternal,lvinternal,extraeqs,(#-1)&/@deltalv,(#-1)&/@integrlv,position/@neededvars],
          FFAlgLinearFitImplemN[gid,inputs,Length[tauvars],lininternal,integrinternal,lvinternal,extraeqs,(#-1)&/@deltalv,(#-1)&/@integrlv,position/@neededvars]
        ]
    ]
];

Options[FFAlgLinearFit] = Join[{"ExtraEquations"->2,"Substitutions"->{}}, Options[FFAlgDenseSolver]];
FFAlgLinearFit[gid_,id_,inputs_List,params_,
                  delta_,integrandin_, tauvarsin_,varsin_,OptionsPattern[]]:=
  FFRegisterAlgorithm[RegisterAlgLinearFit, gid, id, inputs, {params,
                         delta,integrandin, tauvarsin,varsin,
                         OptionValue["ExtraEquations"],OptionValue["Substitutions"],
                         OptionValue["NeededVars"],OptionValue["ApplyFunction"]}];


RegisterAlgSubgraphFit[gid_,inputs_,{subgraphid_,samplevars_,vars_,neededvarsin_,extraeqs_}]:=Module[{neededvars},
  Catch[
    neededvars = If[TrueQ[neededvarsin==Automatic], vars, neededvarsin];
    CheckVariables[samplevars];
    CheckVariables[vars];
    CheckVariables[neededvars];
    FFSubgraphFitImplem[gid,inputs,GetGraphId[subgraphid],Length[vars],Length[samplevars],((Position[vars,#][[1,1]])&/@neededvars)-1,extraeqs]
  ]
];
Options[FFAlgSubgraphFit]={"NeededVars"->Automatic,"ExtraEquations"->2};
FFAlgSubgraphFit[gid_,id_,inputs_List,subgraphid_,samplevars_,coeffs_,OptionsPattern[]]:=FFRegisterAlgorithm[RegisterAlgSubgraphFit,gid,id,inputs,{subgraphid,samplevars,coeffs,OptionValue["NeededVars"],OptionValue["ExtraEquations"]}];


TakeMultiFitElemsToInternal[a___]:>Throw[$Failed];
TakeMultiFitElemsToInternal[a_List->b_List]:=Module[{position},
  If[!SubsetQ[a,Union@@b],Throw[$Failed];];
  position = Association[{}];
  Table[(position[a[[ii]]]=ii-1);,{ii,Length[a]}];
  (position/@#)&/@b
];
RegisterAlgSubgraphMultiFit[gid_,inputs_,{subgraphid_,samplevars_,take_,neededvarsin_,extraeqs_}]:=Module[
  {neededvars},
  Catch[
    neededvars = If[TrueQ[neededvarsin==Automatic], take[[2]], neededvarsin];
    CheckVariables[samplevars];
    CheckVariables[take[[1]]];
    If[!TrueQ[And@@(Table[SubsetQ[take[[2,ii]],neededvars[[ii]]],{ii,Length[take[[2]]]}])],Throw[$Failed]];
    FFSubgraphMultiFitImplem[gid,inputs,GetGraphId[subgraphid],Length[samplevars],
                                TakeMultiFitElemsToInternal[take],
                                Table[((Position[take[[2,ii]],#][[1,1]])&/@neededvars[[ii]])-1,{ii,Length[take[[2]]]}],
                                extraeqs]
  ]
];
Options[FFAlgSubgraphMultiFit]={"NeededVars"->Automatic,"ExtraEquations"->2};
FFAlgSubgraphMultiFit[gid_,id_,inputs_List,subgraphid_,samplevars_,take_,OptionsPattern[]]:=FFRegisterAlgorithm[RegisterAlgSubgraphMultiFit,gid,id,inputs,{subgraphid,samplevars,take,OptionValue["NeededVars"],OptionValue["ExtraEquations"]}];


RegisterAlgSubgraphRec[gid_,inputs_,{subgraphid_,samplevars_,shiftvars_}]:=Module[{neededvars},
  Catch[
    CheckVariables[samplevars];
    FFSubgraphReconstructImplem[gid,inputs,GetGraphId[subgraphid],Length[samplevars],If[TrueQ[shiftvars],1,0]]
  ]
];
Options[FFAlgSubgraphReconstruct]:={"ShiftVars"->True};
FFAlgSubgraphReconstruct[gid_,id_,inputs_List,subgraphid_,samplevars_,OptionsPattern[]]:=FFRegisterAlgorithm[RegisterAlgSubgraphRec,gid,id,inputs,{subgraphid,samplevars,OptionValue["ShiftVars"]}];


FFSubgraphReconstructLearn[gid_,vars_]:=(((Times@@(vars^#))&/@#)&/@#)&/@FFLearn[gid];


FFSubgraphReconstructSol[expr_List,learn_]:=Module[{ppp},
  ppp=PartitionsWithLen[expr,Length[#[[1]]]+Length[#[[2]]]&/@learn];
  Table[((Plus@@#[[1]])/(Plus@@#[[2]]))&[(PartitionsWithLen[ppp[[ii]],{Length[#[[1]]],Length[#[[2]]]}&[learn[[ii]]]]learn[[ii]])],{ii,Length[learn]}]
];


Options[FFLinearFit] := Join[{"IndepVarsOnly"->False},
                                AutoReconstructionOptions[],
                                Options[FFAlgLinearFit]];
FFLinearFit[params_,delta_,integrandin_, tauvarsin_,varsin_,opt:OptionsPattern[]] := Module[
	{res,graph,in,sys,learn,pars,dummy},
	FFNewGraph[graph];
    res = Catch[
      If[!TrueQ[params == {}], CheckVariables[params];];
      FFGraphInputVars[graph,in,params];
      pars = params;
      
      res = FFAlgLinearFit[graph,sys,{in},pars,delta,integrandin,tauvarsin,varsin,
                              Sequence@@FilterRules[{opt}, Options[FFAlgLinearFit]]];
      If[res==$Failed,Throw[$Failed]];
      
      FFGraphOutput[graph,sys];
      FFAlgorithmSetReconstructionOptions[graph,Sequence@@FilterRules[{opt}, Options[FFAlgorithmSetReconstructionOptions]]];
      learn = FFDenseSolverLearn[graph,varsin];
      If[!TrueQ[learn[[0]]==List],Throw[learn]];
      If[TrueQ[OptionValue["IndepVarsOnly"]], Throw["IndepVars"/.learn]];
      
      res = If[TrueQ[params == {}],
               FFReconstructNumeric[graph,Sequence@@FilterRules[{opt}, Options[FFReconstructNumeric]]],
               FFReconstructFunction[graph,params,Sequence@@FilterRules[{opt}, Options[FFReconstructFunction]]]
             ];
       If[!TrueQ[res[[0]]==List],Throw[res]];
       FFDenseSolverSol[res,learn]
    ];
    FFDeleteGraph[graph];
    res
];


FFRatRec[a_List,p_]:=Catch[ToExpression/@FFRatRecImplem[ToString[CheckedInt[#]]&/@a,ToString[CheckedInt[p]]]];
FFRatRec[a_,p_]:=Catch[ToExpression[FFRatRecImplem[{ToString[CheckedInt[a]]},ToString[CheckedInt[p]]][[1]]]];


Options[FFSetLearningOptions]={"PrimeNo"->Automatic,"MaxSingularPoints"->Automatic};
FFSetLearningOptions[graph_,node_,OptionsPattern[]]:=Catch[
FFSetLearningOptionsImplem[
GetGraphId[graph],GetAlgId[graph,node],
Sequence@@(toFFInternalUnsignedFlag[#,OptionValue[#]]&/@Options[FFSetLearningOptions][[;;,1]])]
];


Protect[FFExV];
(* These instructions must be in sync with AnalyticExpression::InstrType in alg_functions.hh *)
FFInstrADD = 0;
FFInstrMUL = 1;
FFInstrNEG = 2;
FFInstrPOW = 3;
FFInstrVAR = 4;
FFInstrNEGPOW = 5;
FFInstrSMALLNUM = 6;
FFInstrMEDNUM = 7;
FFInstrBIGNUM = 8;
FFInstrEND = 9;


RegisterAlgRatExprEval[gid_,inputs_,{params_,funcsin_}]:=Module[
  {fun,FFCompile,numcounter,number,invnumber,error,
   byteinstr,maxsmallint,maxmediumint,bytecodes,funcs},
  
  byteinstr = 2^8;
  maxsmallint = byteinstr-1;
  maxmediumint = 2^62;
  
  funcs = funcsin /. Dispatch[Inner[Rule,params,FFExV/@Range[0,Length[params]-1],List]];
  
  numcounter=0;
  number[i_]:=number[i]=(invnumber[numcounter]=i; numcounter++);
  
  error[expr_]:=(Message[FF::nonratsub,expr/.Dispatch[Inner[Rule,FFExV/@Range[0,Length[params]-1],params,List]]]; Throw[$Failed]);

  FFCompile[expr_]:=error[expr];
  FFCompile[FFExV[j_]]:=(FFCompile[j]; Sow[FFInstrVAR];);
  FFCompile[a_Integer]:=Which[
    0<=a<=maxsmallint,
    Sow[FFInstrSMALLNUM]; Sow[a];,
    0<=-a<=maxsmallint,
    Sow[FFInstrSMALLNUM]; Sow[-a]; Sow[FFInstrNEG];,
    0<=a<=maxmediumint,
    Sow[FFInstrMEDNUM]; (Sow[Length[#]]; Sow/@#;)&@IntegerDigits[a,byteinstr];,
    0<=-a<=maxmediumint,
    Sow[FFInstrMEDNUM]; (Sow[Length[#]]; Sow/@#;)&@IntegerDigits[-a,byteinstr]; Sow[FFInstrNEG];,
    True,
    (FFCompile[number[a]]; Sow[FFInstrBIGNUM];)
  ];
  FFCompile[a_?FFRationalQ]:= (FFCompile[number[a]]; Sow[FFInstrBIGNUM];);
  FFCompile[a_Plus]:=(FFCompile/@(List@@a); FFCompile[Length[a]]; Sow[FFInstrADD]; );
  FFCompile[a_Times]:=(FFCompile/@(List@@a); FFCompile[Length[a]]; Sow[FFInstrMUL];);
  FFCompile[Power[a_,b_Integer]]:=Which[
    TrueQ[b>0],
    (FFCompile[b]; FFCompile[a]; Sow[FFInstrPOW];);,
    TrueQ[b<0],
    (FFCompile[-b]; FFCompile[a]; Sow[FFInstrNEGPOW];);,
    True,
    error[a^b];
  ];
    
  Catch[
    bytecodes = Table[ Reap[FFCompile[fun]; Sow[FFInstrEND];][[2,1]] ,{fun,funcs}];
    FFAlgRatExprEvalImplem[gid,inputs,Length[params],bytecodes,ToString[#,InputForm]&/@(invnumber/@Range[0,numcounter-1])]
  ]  
];

FFAlgRatExprEval[gid_,id_,inputs_List,params_,functions_List]:=FFRegisterAlgorithm[RegisterAlgRatExprEval,gid,id,inputs,{params,functions}];


Options[FFTogether] := Join[{"Parameters"->Automatic},AutoReconstructionOptions[]];
FFTogether[expr_,OptionsPattern[]]:=Module[
  {opt,vars,g,res,in,ex},
  opt = (#[[1]]->OptionValue[#[[1]]])&/@Options[FFTogether];
  vars = OptionValue["Parameters"];
  If[TrueQ[vars==Automatic],
    vars = Variables[expr];
  ];
  If[TrueQ[Length[vars]==0], Return[expr]];
  res = Catch[
    CheckVariables[vars];
    FFNewGraph[g,in,vars];
    If[!TrueQ[FFAlgRatExprEval[g,ex,{in},vars,{expr}]], Throw[$Failed]];
    FFGraphOutput[g,ex];
    If[TrueQ[#[[0]]==List],#[[1]],#]&@FFReconstructFunction[g,vars,Sequence@@FilterRules[{opt}, Options[FFReconstructFunction]]]
  ];
  If[FFGraphQ[g], FFDeleteGraph[g]];
  res
];
FFTogether[expr_List,opt:OptionsPattern[]]:=(FFTogether[#,opt])&/@expr;


fflowlib = $Failed;

FFLoadLib[] := Module[
    {},
    fflowlib = FindLibrary["fflowmlink"];
    If[TrueQ[fflowlib == $Failed],
       Message[FF::nolib];,
       LibraryLoad[fflowlib];
       FFLoadLibObjects[]
    ];
];
FFLoadLibObjects[] := Module[
    {},

    (*FFMulInvImplem = LibraryFunctionLoad[fflowlib, "fflowml_mul_inv", LinkObject, LinkObject];*)

    FFDefaultNThreadsImplem = LibraryFunctionLoad[fflowlib, "fflowml_default_nthreads", LinkObject, LinkObject];

    FFNewGraphImplem = LibraryFunctionLoad[fflowlib, "fflowml_graph_new", LinkObject, LinkObject];
    FFNewDummyGraphImplem = LibraryFunctionLoad[fflowlib, "fflowml_graph_dummy", LinkObject, LinkObject];
    FFDeleteGraphImplem = LibraryFunctionLoad[fflowlib, "fflowml_graph_delete", LinkObject, LinkObject];
    FFDeleteNodeImplem = LibraryFunctionLoad[fflowlib, "fflowml_node_delete", LinkObject, LinkObject];
    FFGraphSetOutputImplem = LibraryFunctionLoad[fflowlib, "fflowml_graph_set_out_node", LinkObject, LinkObject];
    FFSetInputVarsImplem = LibraryFunctionLoad[fflowlib, "fflowml_graph_input_vars", LinkObject, LinkObject];
    FFGraphNParsOutImplem = LibraryFunctionLoad[fflowlib, "fflowml_graph_nparsout", LinkObject, LinkObject];
    FFNodeNParsOutImplem = LibraryFunctionLoad[fflowlib, "fflowml_node_nparsout", LinkObject, LinkObject];
    FFGraphEdgesImplem = LibraryFunctionLoad[fflowlib, "fflowml_graph_edges", LinkObject, LinkObject];
    FFSimpleSubGraphImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_simple_subgraph", LinkObject, LinkObject];
    FFMemoizedSubGraphImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_memoized_subgraph", LinkObject, LinkObject];
    FFSubGraphMapImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_subgraph_map", LinkObject, LinkObject];
    FFSubgraphFitImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_subgraph_fit", LinkObject, LinkObject];
    FFSubgraphMultiFitImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_subgraph_multifit", LinkObject, LinkObject];
    FFSubgraphReconstructImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_subgraph_rec", LinkObject, LinkObject];

    FFRegisterDenseSolverImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_dense_system", LinkObject, LinkObject];
    FFRegisterSparseSolverImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_sparse_system", LinkObject, LinkObject];
    FFRegisterDenseSolverImplemN = LibraryFunctionLoad[fflowlib, "fflowml_alg_num_dense_system", LinkObject, LinkObject];
    FFRegisterSparseSolverImplemN = LibraryFunctionLoad[fflowlib, "fflowml_alg_num_sparse_system", LinkObject, LinkObject];
    FFRegisterNodeDenseSolverImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_node_dense_system", LinkObject, LinkObject];
    FFLearnImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_learn", LinkObject, LinkObject];
    FFSparseSolverMarkAndSweepEqsImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_mark_and_sweep_eqs", LinkObject, LinkObject];
    FFSparseSolverDeleteUnneededEqsImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_delete_unneeded_eqs", LinkObject, LinkObject];
    FFAlgorithmGetInfoImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_get_info", LinkObject, LinkObject];
    FFTotalDegreesImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_degrees", LinkObject, LinkObject];
    FFVarsDegreesImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_vars_degrees", LinkObject, LinkObject];
    FFAllDegreesImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_all_degrees", LinkObject, LinkObject];
    FFSampleImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_sample", LinkObject, LinkObject];
    FFSampleFromPointsImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_sample_from_points", LinkObject, LinkObject];
    FFReconstructFromCurrentEvaluationsImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_reconstruct", LinkObject, LinkObject];
    FFReconstructFromCurrentEvaluationsModImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_reconstruct_mod", LinkObject, LinkObject];
    FFReconstructUnivariateImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_univariate_reconstruct", LinkObject, LinkObject];
    FFReconstructUnivariateModImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_univariate_reconstruct_mod", LinkObject, LinkObject];
    FFReconstructNumericImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_numeric_reconstruct", LinkObject, LinkObject];
    FFAlgMakeNodeMutableImplem = LibraryFunctionLoad[fflowlib, "fflowml_node_set_mutable", LinkObject, LinkObject];
    FFGraphNodesImplem = LibraryFunctionLoad[fflowlib, "fflowml_graph_nodes", LinkObject, LinkObject];
    FFGraphPruneImplem = LibraryFunctionLoad[fflowlib, "fflowml_graph_prune", LinkObject, LinkObject];
    FFSolverResetNeededVarsImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_system_reset_neeed", LinkObject, LinkObject];
    FFSolverOnlyHomogeneousImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_system_only_homogeneous", LinkObject, LinkObject];
    FFSolverSparseOutputImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_system_sparse_output", LinkObject, LinkObject];
    FFIndependentOfImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_independent_of_var", LinkObject, LinkObject];
    FFAlgRatFunEvalImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_ratfun_eval", LinkObject, LinkObject];
    FFAlgRatFunEvalFromCoeffsImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_coeff_ratfun_eval", LinkObject, LinkObject];
    FFAlgRatNumEvalImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_ratnum_eval", LinkObject, LinkObject];
    FFAlgChainImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_chain", LinkObject, LinkObject];
    FFAlgTakeImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_take", LinkObject, LinkObject];
    FFAlgSliceImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_slice", LinkObject, LinkObject];
    FFAlgAddImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_add", LinkObject, LinkObject];
    FFAlgMulImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_mul", LinkObject, LinkObject];
    FFAlgMatMulImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_mat_mul", LinkObject, LinkObject];
    FFAlgSparseMatMulImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_sparse_mat_mul", LinkObject, LinkObject];
    FFAlgNonZeroesImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_nonzero", LinkObject, LinkObject];
    FFAlgTakeAndAddImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_take_and_add", LinkObject, LinkObject];
    FFSolverNIndepEqsImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_n_indep_eqs", LinkObject, LinkObject];
    FFSolverIndepEqsImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_indep_eqs", LinkObject, LinkObject];
    FFAlgJSONSparseSolverImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_json_sparse_system", LinkObject, LinkObject];
    FFAlgJSONRatFunEvalImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_json_ratfun", LinkObject, LinkObject];
    FFAlgLinearFitImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_linear_fit", LinkObject, LinkObject];
    FFAlgLinearFitImplemN = LibraryFunctionLoad[fflowlib, "fflowml_alg_numeric_fit", LinkObject, LinkObject];
    FFAlgLaurentImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_laurent", LinkObject, LinkObject];
    FFDefaultMaxDegImplem = LibraryFunctionLoad[fflowlib, "fflowml_default_maxdeg", LinkObject, LinkObject];
    FFDumpDegreesImplem=LibraryFunctionLoad[fflowlib, "fflowml_dump_degrees", LinkObject, LinkObject];
    FFLoadDegreesImplem=LibraryFunctionLoad[fflowlib, "fflowml_load_degrees", LinkObject, LinkObject];
    FFDumpSamplePointsImplem=LibraryFunctionLoad[fflowlib, "fflowml_dump_sample_points", LinkObject, LinkObject];
    FFDumpEvaluationsImplem=LibraryFunctionLoad[fflowlib, "fflowml_dump_evaluations", LinkObject, LinkObject];
    FFLoadEvaluationsImplem=LibraryFunctionLoad[fflowlib, "fflowml_load_evaluations", LinkObject, LinkObject];
    FFSamplesFileSizeImplem=LibraryFunctionLoad[fflowlib, "fflowml_samples_file_size", LinkObject, LinkObject];
    FFNParsFromDegreesFileImpl=LibraryFunctionLoad[fflowlib, "fflowml_npars_from_degree_info", LinkObject, LinkObject];
    FFGraphEvaluateImplem=LibraryFunctionLoad[fflowlib, "fflowml_graph_evaluate", LinkObject, LinkObject];
    FFGraphEvaluateListImplem=LibraryFunctionLoad[fflowlib, "fflowml_graph_evaluate_list", LinkObject, LinkObject];
    FFPrimeNoImplem=LibraryFunctionLoad[fflowlib, "fflowml_prime_no", LinkObject, LinkObject];
    FFNSamplePointsImplem=LibraryFunctionLoad[fflowlib, "fflowml_alg_count_sample_points", LinkObject, LinkObject];
    FFRatRecImplem=LibraryFunctionLoad[fflowlib, "fflowml_alg_rat_rec", LinkObject, LinkObject];
    FFAlgRatExprEvalImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_ratexpr_eval", LinkObject, LinkObject];
    FFSetLearningOptionsImplem = LibraryFunctionLoad[fflowlib, "fflowml_alg_set_learning_options", LinkObject, LinkObject];
];


FFLoadLib[];


FFAlreadyLoaded = True;


End[] (* "`Private`" *)


EndPackage[] (* "FiniteFlow`" *)
