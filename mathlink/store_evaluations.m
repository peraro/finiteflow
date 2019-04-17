(* ::Package:: *)

(*
 * This shows how to generate sample points and info on total degrees, and save them on disk.
 * Then we can split the evaluations over several machines.
 * In this file I just use FFDeleteGraph to start over, but it should be clear that it means.
 *)


<<FiniteFlow`


(* This code loads an algorithm in a graph *)
LoadAlgorithm[g_]:=(
FFNewGraph[g];
FFGraphInputVars[g,in,{x1,x2}];
FFAlgRatFunEval[g,rf,{in},{x1,x2},{(x1+1)/(x2+2),(x1^4+x2^4+1)/(x1^3-x2^3+27)}];
FFGraphOutput[g,rf];
);


(*
 * This code loads the algorithm and its total degrees, and a list of sample points to be evaluated.
 * The sample points are taken from "samples.fflow".  It starts reading from the sample point "starts"
 * (counting from 0) and it loads a number of points equal to "size".  The evaluations are then
 * stored into "file".
 *)
EvaluateBatch[g_,start_,size_,file_]:=(
  LoadAlgorithm[g];
  FFLoadDegrees[g,"degs.fflow"];
  FFSampleFromPoints[g,"samples.fflow",start,size];  (* or FFAlgorithmSampleFromPoints[g,"samples.fflow",start,size, NTHREADS]; *)
  FFDumpEvaluations[g,file];
);


(* Generate the degrees and a list of sample points. *)


LoadAlgorithm[g]


FFAllDegrees[g]


FFDumpDegrees[g,"degs.fflow"]


(*
 * Note: if we already have some evaluations saved but they weren't enough
 * because MaxPrimes was too low, then we can load them here.  This way the
 * new list of sample points will not include the old ones.
 *) 
(* FFLoadEvaluations[g,{"evals1.fflow","evals2.fflow"}] *)


(* Note: here choose the number of primes you think is needed *)
FFDumpSamplePoints[g,"samples.fflow","MaxPrimes"->1]


FFDeleteGraph[g]


(* read how many points we have generated *)
FFSamplesFileSize["samples.fflow"]


(*
 * We need 26 points.  We divide them into 2 batches of 13 points each.
 * Each batch can be evaluated on a different machine.
 * Then we save the evaluations on disk.
 *)


(* eval batch 1: 13 points starting from point 0 *)


EvaluateBatch[g,0,13,"evals1.fflow"];
FFDeleteGraph[g];


(* eval batch 2: 13 points starting from point 13 *)


EvaluateBatch[g,13,13,"evals2.fflow"];
FFDeleteGraph[g];


(*
 * Here is the reconstruction.  We collect all the evaluation files.
 * We create a "Dummy" graph with 2 inputs and 2 ouputs.
 * Then we load the degrees and the evaluations.  From these we do the reconstruction.
 * If it fails with FFMissingPrimes, it means that more points need to be generated.
 *)


FFNewDummyGraph[d,2,2]


FFLoadDegrees[d,"degs.fflow"]


FFLoadEvaluations[d,{"evals1.fflow","evals2.fflow"}]


(* Note: DO NOT use FFReconstructFunction here *)
FFReconstructFromCurrentEvaluations[d,{x1,x2}]
