# A minimal introduction to the the fflow Python package.

# Load the fflow package and assign to it the "ff" shortcut
import fflow as ff

# NOTE: one may also use "from fflow import *" to import all the
#       definitions of fflow into the public namespace, such that
#       using "ff.*" is not needed to access them. Here we don't do
#       that to make it clear which functions and symbols come from
#       the fflow package.


# We first define a new graph with 3 input variables. This returns
# integer identifiers (id) for the graph and its input node.
graph,innode = ff.NewGraphWithInput(3)

# NOTE: The previous call is a shortcut for:
#         graph = ff.NewGraph()
#         innode = ff.SetGraphInput(graph,3)
# NOTE: The input node always has id=0, but we still save it to a
#       variable for clarity.


# Next we want to define some algorithms in the graph, and combine
# them together in a more complex calculation.  Each algorithm
# corresponds to a node, and gets its inputs from other nodes
# (including the input node defined above).  Most of the commands
# defining algorithms have the form
#
#   newnode = ff.AlgXXX(graph,input1,input2,...,args...)
#
# for algorithms taking a fixed number of input nodes or
#
#   newnode = ff.AlgXXX(graph,[input1,input2,...],args...)
#
# for algorithms taking a variable number of input nodes.  Here XXX is
# the name of the algorithm.  The first argument is the integer ID of
# the graph.  After the first argument, we pass the list of IDs
# corresponding to the input nodes.  Additional arguments and options,
# which are specific to the algorithm XXX, are added at the end.  The
# return value `newnode` is the integer ID of the new node that has
# been created.



# Most problems need some form of analytic input, such as a list of
# rational functions. In this API we have several ways of specifying
# that. The most straightforward one is by letting FiniteFlow parse
# strings with their analytic expressions. The FiniteFlow parser has
# some limitations but it is a viable option (evaluate
# `help(ff.ParseRatFun)` in an interactive shell for more info).
rflist = ff.ParseRatFun(["x1", "x2", "x3"],
                        ["(1+x1)/(1+x2)", "(1+x2)/(1+x3)", "(1+x1)/(1+x3)",
                         "(1+x1^2)/(1+x2^2)"])
# which returns an opaque object that represents a list of rational
# functions. The same kind of object is returned by functional
# reconstruction routines. Additional ways of defining the same list
# are explained later.

# Let us first create a node fun1 which evaluates a list of rational
# functions
fun1 = ff.AlgRatFunEval(graph,innode,rflist)

# Let us create another node fun2 which evaluates another list of
# functions
fun2 = ff.AlgRatFunEval(graph,innode,
                        ff.ParseRatFun("x1,x2,x3".split(","),
                                       "(1+x2^2)/(1+x3^2),\
                                        (1+x1^2)/(1+x3^2),\
                                        (1-x1^2)/(1-x2^2),\
                                        (1-x1^2)/(1-x3^2)".split(",")))

# The nodes fun1 and fun2 return a list of 4 elements each.  We now
# want to interpret this list as the elements of 2*2 matrices and do a
# matrix multiplication between them.
matmul = ff.AlgMatMul(graph,fun1,fun2, 2,2,2)

# We now select this last node as the output node of our dataflow graph.
ff.SetOutputNode(graph,matmul)

# Now the graph can be evaluated over finite fields.  There are a few
# (but not many) cases where you might want to evaluate the graph
# yourself, from Python. In the following we evaluate it with inputs
# [123,345,567] modulo the 0-th prime (whose exact value can be
# checked using `ff.PrimeNo(0)`):
evaluation = ff.EvaluateGraph(graph,[123,345,567],0)

# The output is a list of length 4.  In this case it represents the
# elements of the 2*2 matrix computed via the matrix multiplication
# algorithm.

# In most cases, you will not be interested in evaluating the graph
# yourself, but only in reconstructing the analytic expression of its
# output.  This is easily done with the following command (see also
# `help(ff.ReconstructFunction)` for a list of allowed options):
recfun = ff.ReconstructFunction(graph)

# The reconstructed function `recfun` is an opaque object representing
# a list of functions.  We can check that the reconstruction was
# successful with:
if type(recfun) is ff.RatFunList:
    print("Successful reconstruction!")
else:
    print("Something went terribly wrong here!")
    print("- the reconstruction returned: {}".format(recfun))
    raise ff.ERROR

# We may now convert the list of reconstructed functions to a list of
# strings containing their analytic expression
print("Reconstructed result = ",
      recfun.to_string(["x1", "x2", "x3"]))

# The graph can be deleted once it is no longer needed
ff.DeleteGraph(graph)



# A NOTE about deleting graphs:
#
# FiniteFlow graphs are identified by integers and they must be
# manually deleted using DeleteGraph, when no longer needed.  If a
# graph lives for the whole duration of an application, or close to
# it, deleting it is often unnecessary.  If, however, a graph is only
# needed for a short period of time, failing to delete it will cause a
# memory leak.  For such cases, we offer the utilities GraphContext
# and GraphContextWithInput which create a context in which the graph
# lives and automatically delete the graph when we exit the block
# which defines it (even in case of an error or exception).
#
# Examples of usage are:
with ff.GraphContext() as graph:
    ratnum = ff.AlgRatNumEval(graph, ["3","2","1"])
    ff.SetOutputNode(graph, ratnum)
    print("Graph evaluation yields: ", ff.EvaluateGraph(graph,[],0))
# and
with ff.GraphContextWithInput(2) as (graph,inputnode):
    ratfun = ff.AlgRatFunEval(graph, inputnode,
                              ff.ParseRatFun(["x","y"],["x+y","x^2+1"]))
    ff.SetOutputNode(graph, ratfun)
    print("Graph evaluation yields: ", ff.EvaluateGraph(graph,[3,2],0))
# In both cases the graph is deleted as soon as the "with ..." block
# is exited.


# FINAL NOTES:
#
# We conclude by listing additional methods for defining lists of
# rational functions. One of them consists in defining them via their
# monomial data, namely a list of coefficients and monomial exponents
# for numerator and denominator.
mondata = [
    # (1+3/2 x1^2 x2)/(1-x2)
    ([('1', (0, 0, 0)), ('3/2', (2, 1, 0))],
     [('1', (0, 0, 0)), ('-1', (0, 1, 0))]),
    # (-1/14+x1 x2^2 x3^3)/(-7+x3)
    ([('-1/14', (0, 0, 0)), ('1', (1, 2, 3))],
     [('-7', (0, 0, 0)), ('1', (0, 0, 1))]),
]
rflist = ff.NewRatFunList(3, mondata)
# which should return the same as
rflist2 = ff.ParseRatFun(["x1","x2","x3"],
                         ["(1+3/2 x1^2 x2)/(1-x2)",
                          "(-1/14+x1 x2^2 x3^3)/(-7+x3)"])

# The monomial data for any RatFunList can be retrieved using the
# monomials() method. We can check that it is identical for the two
# functions we just defined:
if not (rflist.monomials() == rflist2.monomials() == mondata):
    print("Monomial data is not what we expected!")
    raise ff.ERROR

# Finally, we can JSON-serialize such a list of functions using:
#
#   ff.RatFunToJSON(rflist, "filename.json")
#
# and define a node that evaluates those using
#
#   ff.AlgJSONRatFunEval(graph,innode,"filename.json")
#
# which might be faster than other options, since it avoids passing
# monomial data from/to Python.
