#include <fflow/capi.h>
#include <stdio.h>

// Note: on failure we are leaking memory, but it's not worth fixing
// for a simple test like this

FFStatus test_jsonratfun()
{
  // Define graph from JSON-serialized list of rational functions

  const char * input_file = "tests_json/test_ratfuns.json";
  const unsigned nvars = 3;

  FFNode in = 0;
  FFGraph graph = ffNewGraphWithInput(nvars, &in);
  if (ffIsError(graph)) {
    fprintf(stderr, "Failed to create graph with %u input variables.\n", in);
    return FF_ERROR;
  }

  FFNode rf = ffAlgJSONRatFunEval(graph, in, input_file);
  if (ffIsError(rf)) {
    fprintf(stderr, "Failed rf node from JSON file %s.\n", input_file);
    return FF_ERROR;
  }

  if (ffSetOutputNode(graph, rf) != FF_SUCCESS) {
    fprintf(stderr, "Failed to set output node.");
    return FF_ERROR;
  }

  FFRatFunList * results;
  FFStatus ret = ffReconstructFunction(graph,
                                       (FFRecOptions){.max_primes = 10},
                                       &results);

  if (ret != FF_SUCCESS) {
    fprintf(stderr, "Failed to reconstruct function: ");
    switch(ret) {
    case FF_MISSING_POINTS:
      fprintf(stderr, "missing points.\n");
      break;
    case FF_MISSING_PRIMES:
      fprintf(stderr, "missing primes.\n");
      break;
    default:
      fprintf(stderr, "unknown error.\n");
    }
    return FF_ERROR;
  }

  // Evaluate graph at a random-like point
  const FFUInt x[3] = {3057312585776011302ULL,
                       3795153781312484964ULL,
                       3415194000889226426ULL};
  const unsigned prime_no = 10;
  FFUInt * output1 = ffEvaluateGraph(graph, x, prime_no);
  if (!output1) {
    fprintf(stderr, "Failed to evaluate graph for check.\n");
    return FF_ERROR;
  }

  // Compare with numerical evaluation of result

  FFUInt * output2 = ffEvaluateRatFunList(results, x, prime_no);
  if (!output2) {
    fprintf(stderr, "Failed to evaluate graph for check.\n");
    return FF_ERROR;
  }

  const unsigned nparsout = ffGraphNParsOut(graph);
  ret = FF_SUCCESS;
  for (unsigned j=0; j<nparsout; ++j)
    if (output1[j] != output2[j]) {
      fprintf(stderr, "Failed comparison at %u-th position.\n", j);
      ret = FF_ERROR;
    }

  if (ret == FF_SUCCESS)
    printf("Successful reconstruction!\n");

  ffFreeMemoryU64(output1);
  ffFreeMemoryU64(output2);
  ffFreeRatFun(results);

  return ret;
}


int main(void)
{
  return test_jsonratfun();
}
