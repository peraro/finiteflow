# In this example we cover algorithms which have a learning phase.  In
# particular, we focus on linear systems.
#
# This example is similar to the one in `1_basic_usage.py`, but this
# time we want to compute the inverse of a matrix, and then we use it
# in a matrix multiplication.  We do not want to reconstruct the
# inverse analytically, but only the result of the final matrix
# multiplication.


import fflow as ff

# We create a new graph with three input variables
xs = ["x1", "x2", "x3"]
graph,innode = ff.NewGraphWithInput(len(xs))

# Define two 2x2 analytic matrices as a list of their 4 entries in
# row-major order
mat1 = "(1+x2^2)/(1+x3^2), (1+x1^2)/(1+x3^2), (1-x1^2)/(1-x2^2), (1-x1^2)/(1-x3^2)".split(",")
mat2 = "(1+x1)/(1+x2), (1+x2)/(1+x3), (1+x1)/(1+x3), (1+x1^2)/(1+x2^2)".split(",")

# We want to compute "mat1 . inverse(mat2)"

# First, we define a node computing the entries of mat1:
m1 = ff.AlgRatFunEval(graph,innode,ff.ParseRatFun(xs,mat1))

# As for the inverse of mat2, we can define it, using the Gauss-Jordan
# method, as the solution of the following linear system:
#
#   \sum_j mat2(i,j)*z_j - t_i == 0  for i=1,2
#
# with respect to the unknowns (sorted by weight) {z_1, z_2, t_1,
# t_2}.  Notice that the ordering is important.
#
# We define the system of equations as a sparse system from its matrix
# of coefficients.  The system
#
#    A.x == b,
#
# where A is a matrix, x a vector of unknowns, and b a known vector is
# defined by the matrix
#
#   (A|b).
#
# For each row we need to specify the list of coefficients and the
# list of nonzero columns (which correspond to the unknowns
# multiplying the coefficients):
Aij = [
    # non-zero coefficients of first and second row of:
    #  mat2.(z1,z2) - (t1,t2) == 0
    mat2[0], mat2[1], "-1",
    mat2[2], mat2[3], "-1"
]
cols = [
    # list of unknowns in each equation, as indexes in the list
    # [z1,z2,t1,t2]
    [0,1,2],
    [0,1,3]
]
ccs = ff.ParseIdxRatFun(xs, Aij)
n_unknowns = 4
invmat2 = ff.AlgAnalyticSparseLSolve(graph, innode, n_unknowns, cols, ccs)

# We don't want the constant term in the solution (which is zero), so
# we exclude it with
ff.LSolveOnlyHomogeneous(graph,invmat2)

# A linear solver needs to complete a learning phase, where it learns
# the dependent variables, independent variables, etc...  Before this
# is completed, the node cannot be evaluated.
ff.SetOutputNode(graph,invmat2)
ff.Learn(graph)

# After the learning stage, we can access info about dependent and
# independent variables.  From thes we deduce that the matrix mat is
# invertible.
if ff.LSolveDepVars(graph,invmat2) != [0,1]:
    print("Unexpected list of dependent variables!")
if ff.LSolveIndepVars(graph,invmat2) != [2,3]:
    print("Unexpected list of independent variables!")
    raise ff.FFlowError()

# NOTE: LSolveDepVars(...) are not guaranteed to be sorted, so in
#       general we may need to reorder the output before moving to the
#       next step, which we may do using an ff.AlgTake node.

# Thes next two calls don't do anything useful (nor harmful) here, but
# they are generally recommended to optimize sparse linear systems.
ff.LSolveMarkAndSweepEqs(graph, invmat2)
ff.LSolveDeleteUnneededEqs(graph, invmat2)

# We can (optionally) check that invmat2 is indeed computing the
# inverse of mat2, by creating a node for mat2 and one which
# multiplies it with its inverse.  This part is optional.

check_inverse = True # set to False to skip the check

if check_inverse:

    # define mat2 matrix
    m2 = ff.AlgRatFunEval(graph, innode, ff.ParseRatFun(xs,mat2))

    # compute inverse(mat2) . mat2
    checknode = ff.AlgMatMul(graph,invmat2,m2, 2,2,2)

    # evaluate the result (numerically) at some random-like point to
    # check it is the (flattened out) identity matrix
    ff.SetOutputNode(graph,checknode)
    checkeval = ff.EvaluateGraph(graph,[12345678,3456789,456789678],0)
    idmat = [1,0,
             0,1]
    if checkeval == idmat:
        print("Successfully checked: `inverse(mat2) . mat2 == identity`")
    else:
        print("Check of inverse matrix failed!")
        raise ff.FFlowError()

# We now perform the matrix multiplication we said we wanted,
# i.e. mat1.invmat2...
matmul = ff.AlgMatMul(graph,m1,invmat2, 2,2,2)

# ...and reconstruct its output
ff.SetOutputNode(graph,matmul)
rec = ff.ReconstructFunction(graph)
ff.DeleteGraph(graph)

# This prints a list of strings with the expressions of the entries of
# the reconstructed matrix (flattened out)
print("Reconstructed : ", rec.to_string(xs))
