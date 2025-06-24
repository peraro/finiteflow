'''
Adapted from ../test_sparse_solver.wl
'''

from fflow import *

cols = [[11, 17, 23, 29], [10, 16, 22, 28], [17, 23, 35], [5, 9, 11, 15, 16, 17, 21, 22, 27, 34], [4, 10, 16], [15, 21, 33], [3, 9, 15], [11, 29, 35], [5, 8, 10, 11, 14, 17, 20, 26, 28, 34], [4, 10, 16], [5, 9, 11, 14, 17, 20, 27, 32, 33], [2, 3, 4, 8, 9, 10, 14, 15, 16], [3, 9, 15], [8, 26, 32], [2, 8, 14], [2, 8, 14], [11, 29, 35], [5, 7, 10, 13, 19, 23, 25, 28, 29, 34], [4, 22, 28], [5, 9, 11, 13, 19, 23, 27, 31, 33, 35], [1, 3, 4, 5, 7, 10, 13, 21, 22, 27, 34], [4], [3, 9, 21, 33], [3], [7, 8, 11, 25, 26, 29, 31, 32, 35], [1, 2, 5, 7, 10, 13, 20, 26, 28, 34], [4], [1, 2, 5, 7, 8, 9, 13, 20, 27, 32, 33], [2, 3, 4], [3], [8, 26, 32], [2], [2], [7, 25, 31], [1, 19, 25], [1, 7, 19, 31], [1], [7, 25, 31], [1], [1], [17, 23, 35], [5, 6, 12, 16, 18, 22, 23, 24, 29, 34], [4, 22, 28], [12, 15, 17, 18, 21, 23, 30, 33, 35], [0, 3, 5, 6, 12, 16, 21, 22, 27, 34], [4], [15, 21, 33], [3], [5, 6, 14, 17, 20, 24, 29, 30, 32, 35], [0, 2, 4, 5, 6, 12, 16, 20, 26, 28, 34], [4], [0, 3, 5, 6, 12, 14, 15, 20, 27, 32, 33], [2, 3, 4], [3], [2, 14, 26, 32], [2], [2], [5, 6, 13, 19, 23, 24, 29, 30, 31], [0, 1, 4, 18, 19, 22, 24, 25, 28], [0, 3, 5, 6, 13, 18, 19, 21, 27, 30, 31], [0, 1, 4], [3], [1, 2, 5, 6, 13, 20, 24, 25, 26, 30, 31], [0, 1, 4], [0, 1, 2, 3], [2], [1, 19, 25], [1], [1], [12, 18, 30], [0, 18, 24], [12, 18, 30], [0], [0, 12, 24, 30], [0], [0], [0, 18, 24], [0], [0]]
coeffs = ["y", "y", "y", "y", "y", "y", "y", "y", "y", "y", "y", "1", "y", "y", "y", "y", "y", "y", "y", "y", "y", "1", "y", "y", "y", "y", "y", "1", "y", "y", "y", "y", "y", "1", "y", "y", "y", "y", "y", "y", "y", "y", "y", "1", "y", "y", "1", "y", "y", "y", "y", "y", "y", "y", "y", "1", "1", "1", "y", "y", "y", "y", "y", "y", "1", "y", "y", "y", "y", "y", "1", "y", "y", "1", "y", "y", "y", "y", "y", "1", "y", "y", "y", "y", "y", "y", "y", "y", "y", "1", "y", "y", "1", "y", "y", "y", "y", "y", "y", "y", "y", "y", "1", "1", "1", "1", "y", "y", "y", "y", "y", "y", "y", "1", "1", "y", "y", "y", "1", "y", "y", "y", "y", "y", "y", "y", "y", "y", "1", "1", "1", "y", "y", "y", "y", "y", "y", "y", "1", "1", "1", "1", "y", "y", "y", "y", "y", "y", "y", "y", "1", "1", "1", "1", "y", "y", "y", "1", "1", "y", "y", "y", "1", "y", "y", "1", "y", "y", "y", "1", "y", "y", "y", "1", "1", "y", "y", "y", "1", "y", "y", "y", "y", "y", "y", "y", "y", "y", "1", "y", "y", "y", "y", "y", "y", "y", "y", "y", "y", "y", "1", "1", "1", "y", "y", "y", "y", "y", "y", "y", "1", "y", "y", "y", "1", "1", "y", "y", "y", "y", "y", "y", "y", "y", "y", "1", "1", "1", "1", "y", "y", "y", "y", "y", "y", "y", "1", "1", "1", "1", "y", "y", "y", "y", "y", "y", "y", "y", "1", "1", "1", "1", "1", "y", "y", "y", "1", "1", "1", "y", "y", "y", "y", "y", "y", "y", "y", "1", "1", "1", "y", "y", "y", "y", "y", "y", "1", "1", "1", "y", "y", "y", "y", "y", "y", "y", "y", "1", "1", "1", "1", "1", "1", "1", "y", "y", "y", "y", "y", "y", "y", "y", "1", "1", "1", "1", "1", "1", "1", "1", "1", "y", "y", "1", "1", "y", "y", "y", "1", "y", "y", "y", "y", "y", "1", "1", "y", "y", "y", "1", "1", "1", "y", "y", "1", "1"]
depvars = [15, 35, 10, 8, 14, 9, 4, 22, 3, 16, 2, 1, 19, 25, 28, 0, 18, 24, 26, 21, 17, 12, 11, 7, 5, 20, 13, 23, 6, 27, 29, 30]


def testSolve():
    print("Test solver...")
    g,inp = NewGraphWithInput(1)
    ls = AlgAnalyticSparseLSolve(g,inp,36, cols, ParseIdxRatFun(["y"],coeffs))
    SetOutputNode(g,ls)
    Learn(g)
    if LSolveDepVars(g,ls) != depvars:
        print("- Dep vars = ",LSolveDepVars(g,ls))
        print("- Test failed: something wrong with the system")
        exit(1)
    print("- Test passed")


def testMaxCol(maxcol,backsubst,keepfullout):
    print("Test solver max col with inputs: {}, {}, {}"
          .format(maxcol,backsubst,keepfullout))

    ccs = ["1", "2", "3", "a",
           "7", "3", "-4", "b",
           "1", "3", "13", "c"]
    cols = [(0,1,2,3) for _ in range(3)]
    params = ["a","b","c"]

    g,inp = NewGraphWithInput(len(params))
    ls = AlgAnalyticSparseLSolve(g,inp,3,cols,ParseIdxRatFun(params,ccs),
                                 needed_vars=(1,2))

    LSolveSparseOutputWithMaxCol(g,ls,maxcol,backsubst,keepfullout)

    SetOutputNode(g,ls)
    Learn(g)
    LSolveMarkAndSweepEqs(g,ls)
    LSolveDeleteUnneededEqs(g,ls)

    rec = ReconstructFunction(g)
    print("- Reconstructed: " + str(rec.to_string(params)))

if __name__ == '__main__':
    testSolve()
    testMaxCol(1,False,True)
    testMaxCol(2,False,True)
    testMaxCol(1,True,True)
    testMaxCol(2,True,True)
    testMaxCol(1,True,False)
    testMaxCol(2,True,False)
    testMaxCol(1,False,False)
    testMaxCol(2,False,False)
