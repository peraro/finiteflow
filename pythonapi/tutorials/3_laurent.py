# In this tutorial we introduce subgraph nodes.  In particular, we
# focus on the Laurent series expansion algorithm.

import fflow as ff

# We define a first graph g1 depending on three variables.  To make
# things simple, we just define a list of rational functions, but it
# could be any other arbitrarily complicated graph, whose output may
# not be known analytically.

xs = ["x1", "x2", "x3"]
g1,innode1 = ff.NewGraphWithInput(len(xs))

toexpand = ff.ParseRatFun(xs, "(1+x1 x2^2)/(x1+x1 x3^2), (1+x1^2)/(x1+x1^2 x3^2), (x1-x1^3)/(1-x1 x2^2), (1-x1^2)/(x1-x3^2)".split(","))

functions = ff.AlgRatFunEval(g1,innode1,toexpand)
ff.SetOutputNode(g1, functions)

# We now want to define a graph g2 which computes the Laurent
# expansion of g1 with respect to x1 (which must be the first input
# variable of g1) up to (and including) x^order terms, where order can
# be any positive, zero, or negative integer, say

order=3

# NOTE: We can also specify a different order for each function to be
#       Laurent-expanded by making `order` a list instead

# The coefficients of the Laurent expansion will be rational functions
# of x2 and x3. We create the graph g2 that evaluates them.

newxs = ["x2", "x3"]
g2,innode2 = ff.NewGraphWithInput(len(newxs))

laurent = ff.AlgLaurent(g2,innode2,g1,order)
ff.SetOutputNode(g2,laurent)

# Laurent expansion algorithms need a learning phase, where the
# starting order of the expansion and other information is obtained.

ff.Learn(g2)
starting_orders = ff.LaurentMinOrders(g2,laurent)

# We can easily generate, for each function in the output of the graph
# g1, a list of pairs (o,j), meaning that the coefficient of O(x1^o)
# in the Laurent expansion is the j-th element of the output of graph
# g2
order_idx_pairs = ff.LaurentOutput(g2,laurent,range(ff.GraphNParsOut(g2)))


# The coefficients of the Laurent expansion can be computed
# numerically and reconstructed.
recratfun = ff.ReconstructFunction(g2)

# Similarly as before, we may list, for each function we have
# expanded, a list of pairs (o,c), meaning that c is the coefficient
# (converted to a string) of O(x1^o) of the expansion
order_coeff_pairs = ff.LaurentOutput(g2,laurent,recratfun.to_string(newxs))

# Notice that, as mentioned, the coefficients of the Laurent expansion
# computed by g2 can be obtained without knowing the analytic
# expression of the output of g1.

# Cleaning up
ff.DeleteGraph(g2)
ff.DeleteGraph(g1)

# "Nicer" printing of results
print("\n")
for f in order_coeff_pairs:
    print(" + ".join("( {} )*x1^{}".format(c,e) for e,c in f) +
          " + O(x1^{})\n".format(order+1))
