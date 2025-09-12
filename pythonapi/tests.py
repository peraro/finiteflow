from fflow import *


def testRatFun():
    print("Testing rational functions...")

    import pathlib
    thisfile = pathlib.Path(__file__).resolve()

    input_file = str(thisfile.parent.parent.resolve() / "tests_json" / "test_ratfuns.json")

    # build graph
    nvars = 3
    graph,innode = NewGraphWithInput(nvars)
    rf = AlgJSONRatFunEval(graph, innode, input_file)
    SetOutputNode(graph, rf)

    # reconstruct
    rec = ReconstructFunction(graph, max_primes = 10)

    # testing at random-like point
    x = [3057312585776011302, 3795153781312484964, 3415194000889226426]
    prime_no = 10
    output1 = EvaluateGraph(graph, x, prime_no)

    # compare with numerical evaluation of result
    output2 = EvaluateRatFunList(rec, x, prime_no)

    if output1 != output2:
        print("- Test 1/3 failed :(")
        exit(1)
    print("- Test 1/3 passed!!!")

    # build another graph evaluating the reconstructed result
    graph2,innode2 = NewGraphWithInput(nvars)
    rf2 = AlgRatFunEval(graph2, innode2, rec)
    SetOutputNode(graph2, rf2)
    output3 = EvaluateGraph(graph2, x, prime_no)

    # compare again
    if output1 != output3:
        print("- Test 2/3 failed :(")
        exit(1)
    print("- Test 2/3 passed!!!")


    degdata = ParallelReconstructDegreeData(graph)
    rec2 = ReconstructFunctionWithDegrees(graph,degdata)
    if rec2.monomials() != rec.monomials():
        print("- Test 3/3 failed :(")
        exit(1)
    DeleteGraph(graph)
    DeleteGraph(graph2)
    print("- Test 3/3 passed!!!")


def testParsing():
    print("Testing parsing of rational functions...")
    funcs = "(1 + x1 + x1^60 + x1^61 + x2 + x1^60*x2 + x2^60 + x1*x2^60 + x2^61 + x3^60 + x1*x3^60 + x2*x3^60)/(1 + 256*x1 + 18446744073709551616*x2),  256 + x1 + x2 + x3, (1 + x1 + x1^60 + x1^61 + x2 + x1^60*x2 + x2^60 + x1*x2^60 + x2^61 + x3^60 + x1*x3^60 + x2*x3^60)/  (2 + 65536*x1 + 18446744073709551616*x2), 65536 + x1 + x2 + x3,  (1 + x1 + x1^60 + x1^61 + x2 + x1^60*x2 + x2^60 + x1*x2^60 + x2^61 + x3^60 + x1*x3^60 + x2*x3^60)/(3 + 16777216*x1 + 18446744073709551616*x2),  16777216 + x1 + x2 + x3, (1 + x1 + x1^60 + x1^61 + x2 + x1^60*x2 + x2^60 + x1*x2^60 + x2^61 + x3^60 + x1*x3^60 + x2*x3^60)/  (4 + 4294967296*x1 + 18446744073709551616*x2), 4294967296 + x1 + x2 + x3,  (1 + x1 + x1^60 + x1^61 + x2 + x1^60*x2 + x2^60 + x1*x2^60 + x2^61 + x3^60 + x1*x3^60 + x2*x3^60)/(5 + 1099511627776*x1 + 18446744073709551616*x2), 1099511627776 + x1 + x2 + x3,  (1 + x1 + x1^60 + x1^61 + x2 + x1^60*x2 + x2^60 + x1*x2^60 + x2^61 + x3^60 + x1*x3^60 + x2*x3^60)/(6 + 281474976710656*x1 + 18446744073709551616*x2), 281474976710656 + x1 + x2 + x3,  (1 + x1 + x1^60 + x1^61 + x2 + x1^60*x2 + x2^60 + x1*x2^60 + x2^61 + x3^60 + x1*x3^60 + x2*x3^60)/(7 + 72057594037927936*x1 + 18446744073709551616*x2), 72057594037927936 + x1 + x2 + x3,  (1 + x1 + x1^60 + x1^61 + x2 + x1^60*x2 + x2^60 + x1*x2^60 + x2^61 + x3^60 + x1*x3^60 + x2*x3^60)/(8 + 18446744073709551616*x1 + 18446744073709551616*x2), 18446744073709551616 + x1 + x2 + x3,  (1 + x1 + x1^60 + x1^61 + x2 + x1^60*x2 + x2^60 + x1*x2^60 + x2^61 + x3^60 + x1*x3^60 + x2*x3^60)/(9 + 4722366482869645213696*x1 + 18446744073709551616*x2), 4722366482869645213696 + x1 + x2 + x3,  (1 + x1 + x1^60 + x1^61 + x2 + x1^60*x2 + x2^60 + x1*x2^60 + x2^61 + x3^60 + x1*x3^60 + x2*x3^60)/(10 + 1208925819614629174706176*x1 + 18446744073709551616*x2), 1208925819614629174706176 + x1 + x2 + x3"
    funcs = funcs.split(",")
    check = [8016150761108817393, 1044288331122947657, 4086451579588194672, 1044288331123012937, 9124275284856115687, 1044288331139724617, 4114689361543759826, 1044288335417914697, 8146000013182627095, \
1044289430634575177, 5570016077230430187, 1044569806099658057, 4028239444661437179, 1116345925160875337, 1559219309039818263, 1044288331122948435, 2709931136303553563, 1044288331123212105, 2521090694700537326, 1044288331190711625]

    rf = ParseRatFun(["x1", "x2", "x3"], funcs)
    x = [3057312585776011302, 3795153781312484964, 3415194000889226426]
    prime_no = 10
    res = EvaluateRatFunList(rf,x,prime_no)

    if res != check:
        print("- Test failed :(")
        exit(1)
    print("- Test passed!!!")


def testTutorial2():
    print("Test (variant of) second example of Mathematica tutorial...")
    import pathlib, json
    thisfile = pathlib.Path(__file__).resolve()

    # Flattened mat1 and mat2 matrices
    mat1 = ["(1 + x2^2)/(1 + x3^2)",
            "(1 + x1^2)/(1 + x3^2)",
            "(1 - x1^2)/(1 - x2^2)",
            "(1 - x1^2)/(1 - x3^2)"]
    mat2 = ["(1 + x1)/(1 + x2)",
            "(1 + x2)/(1 + x3)",
            "(1 + x1)/(1 + x3)",
            " (1 + x1^2)/  (1 + x2^2)"]
    params = ["x1", "x2", "x3"]

    # The system for the inverse of mat2 is encoded in JSON format
    json_file = str(thisfile.parent.parent.resolve() / "tests_json" / "tinvsys.json")
    json_eqsfile = str(thisfile.parent.parent.resolve() / "tests_json" / "tinveqs.json")
    with open(json_file, 'w') as f:
        eqs_json = [2,4,3,4,[0,1,2,3],1,[json_eqsfile]]
        json.dump(eqs_json, f)


    mygraph, myinput = NewGraphWithInput(3)

    m1rf = ParseRatFun(params, mat1)
    m1 = AlgRatFunEval(mygraph, myinput, m1rf)

    invmat2 = AlgJSONSparseLSolve(mygraph, myinput, json_file)
    LSolveOnlyHomogeneous(mygraph, invmat2)
    SetOutputNode(mygraph, invmat2)
    Learn(mygraph)
    depvars = LSolveDepVars(mygraph, invmat2)
    indepvars = LSolveIndepVars(mygraph, invmat2)
    if (depvars != [0,1] or indepvars != [2,3]):
        print("- Test failed: could not invert the matrix")
        exit(1)

    # checking the inverse works
    m2rf = ParseRatFun(params, mat2)
    m2 = AlgRatFunEval(mygraph, myinput, m2rf)
    check = AlgMatMul(mygraph, invmat2, m2, 2, 2, 2)
    SetOutputNode(mygraph, check)
    if (EvaluateGraph(mygraph, [123,456,789], 0) != [1,0,0,1]):
        print("- Test failed: wrong inverse of matrix")
        exit(1)

    # computing mat1.inverse(mat2)
    matmul = AlgMatMul(mygraph, m1, invmat2, 2, 2, 2)
    SetOutputNode(mygraph, matmul)
    rec = ReconstructFunction(mygraph)

    # checking...
    x = [3057312585776011302, 3795153781312484964, 3415194000889226426]
    prime_no = 10
    numerical_check = [5997495238819826275, 7935486666042592923,
                       4481828051701803585, 5853523369524799279]
    analytic_check = ParseRatFun(params, testTutorial2AnalyticCheck())

    if EvaluateRatFunList(analytic_check, x, prime_no) != numerical_check:
        print("- Test failed: error in parsing analytic result")
        exit(1)

    if EvaluateRatFunList(rec, x, prime_no) != numerical_check:
        print("- Test failed: reconstructed wrong result")
        exit(1)

    DeleteGraph(mygraph)
    print("- Test passed!")



def testTutorial2AnalyticCheck():
    return '''
    (-x1 - x1^3 - x1*x2 - x1^3*x2 - x1*x2^2 - x1^3*x2^2 - x1*x2^3 -
   x1^3*x2^3 + x3 - x1*x3 + x1^2*x3 - x1^3*x3 + x2*x3 - x1*x2*x3 +
   x1^2*x2*x3 - x1^3*x2*x3 + x2^2*x3 - x1*x2^2*x3 + x1^2*x2^2*x3 -
   x1^3*x2^2*x3 + x2^3*x3 - x1*x2^3*x3 + x1^2*x2^3*x3 - x1^3*x2^3*x3 +
   x3^2 + x1^2*x3^2 + x2*x3^2 + x1^2*x2*x3^2 + x2^2*x3^2 + x1^2*x2^2*x3^2 +
   x2^3*x3^2 + x1^2*x2^3*x3^2)/(x1^2 + x1^3 - 2*x2 - 2*x1*x2 - 2*x2^2 -
   2*x1*x2^2 - 2*x2^3 - 2*x1*x2^3 - x2^4 - x1*x2^4 + 2*x3 + 2*x1*x3 +
   2*x1^2*x3 + 2*x1^3*x3 + x3^2 + x1*x3^2 + 2*x1^2*x3^2 + 2*x1^3*x3^2 -
   2*x2*x3^2 - 2*x1*x2*x3^2 - 2*x2^2*x3^2 - 2*x1*x2^2*x3^2 - 2*x2^3*x3^2 -
   2*x1*x2^3*x3^2 - x2^4*x3^2 - x1*x2^4*x3^2 + 2*x3^3 + 2*x1*x3^3 +
   2*x1^2*x3^3 + 2*x1^3*x3^3 + x3^4 + x1*x3^4 + x1^2*x3^4 + x1^3*x3^4),
 (x1 + x1^2 + x1^3 - 2*x2 - 2*x2^2 + x1*x2^2 + x1^2*x2^2 + x1^3*x2^2 -
   4*x2^3 - 3*x2^4 - 2*x2^5 - x2^6 + x3 + 2*x1*x3 + 2*x1^2*x3 + 2*x1^3*x3 -
   2*x2*x3 - x2^2*x3 + 2*x1*x2^2*x3 + 2*x1^2*x2^2*x3 + 2*x1^3*x2^2*x3 -
   4*x2^3*x3 - 3*x2^4*x3 - 2*x2^5*x3 - x2^6*x3 + x3^2 + x1*x3^2 +
   x1^2*x3^2 + x1^3*x3^2 + x2^2*x3^2 + x1*x2^2*x3^2 + x1^2*x2^2*x3^2 +
   x1^3*x2^2*x3^2)/(x1^2 + x1^3 - 2*x2 - 2*x1*x2 - 2*x2^2 - 2*x1*x2^2 -
   2*x2^3 - 2*x1*x2^3 - x2^4 - x1*x2^4 + 2*x3 + 2*x1*x3 + 2*x1^2*x3 +
   2*x1^3*x3 + x3^2 + x1*x3^2 + 2*x1^2*x3^2 + 2*x1^3*x3^2 - 2*x2*x3^2 -
   2*x1*x2*x3^2 - 2*x2^2*x3^2 - 2*x1*x2^2*x3^2 - 2*x2^3*x3^2 -
   2*x1*x2^3*x3^2 - x2^4*x3^2 - x1*x2^4*x3^2 + 2*x3^3 + 2*x1*x3^3 +
   2*x1^2*x3^3 + 2*x1^3*x3^3 + x3^4 + x1*x3^4 + x1^2*x3^4 + x1^3*x3^4),
 (-x1 + 2*x1^2 - x1^3 + x2^4 - x1^2*x2^4 + x3 - x1*x3 + x1^2*x3 - x1^3*x3 -
   x3^2 + x1*x3^2 - x1^2*x3^2 + x1^3*x3^2 - x3^3 + x1*x3^3 - x1^2*x3^3 +
   x1^3*x3^3)/(x1^2 - 2*x2 - x1^2*x2 + x2^4 + x2^5 + 2*x3 + x1^2*x3 -
   x1^2*x2*x3 - x2^4*x3 - x2^5*x3 - x3^2 - x1^2*x3^2 + x2*x3^2 +
   x1^2*x2*x3^2 - x3^3 - x1^2*x3^3 + x2*x3^3 + x1^2*x2*x3^3),
 (x1 - x1^2 - 2*x2 + x1*x2 + x1^2*x2 + x1*x2^2 - x1^2*x2^2 - 2*x2^3 +
   x1*x2^3 + x1^2*x2^3 + x3 - x1^2*x3 - x2*x3 + x1^2*x2*x3 + x2^2*x3 -
   x1^2*x2^2*x3 - x2^3*x3 + x1^2*x2^3*x3 + x3^2 - x1*x3^2 + x2*x3^2 -
   x1*x2*x3^2 + x2^2*x3^2 - x1*x2^2*x3^2 + x2^3*x3^2 - x1*x2^3*x3^2)/
  (x1^2 - 2*x2 - x1^2*x2 + x2^4 + x2^5 + 2*x3 + x1^2*x3 - x1^2*x2*x3 -
   x2^4*x3 - x2^5*x3 - x3^2 - x1^2*x3^2 + x2*x3^2 + x1^2*x2*x3^2 - x3^3 -
   x1^2*x3^3 + x2*x3^3 + x1^2*x2*x3^3)
    '''.split(",")



def testBasicRatFunInterface():
    print("Test rational function interfaces")

    # From coefficiensts and exponents
    ratfunlist1 = [
        ([("1/2", [1,0,2]), ("2/3",[0,2,3])], # numerator
         [("1", [0,0,0]), ("-5/7",[1,0,3])]   # denominator
         ),
        ([("-2", [0,0,0]), ("8",[21,22,23]), ("-729/92",[41,42,43])], # num.
         [("1", [0,0,0])] # den.
         )
    ]
    # from analytic expression
    ratfunlist2 = [
        "(1/2 x1 x3^2 + 2/3 x2^2 x3^3)/(1 - 5/7 x1 x3^3)",
        "-2 + 8 x1^21 x2^22 x3^23 - 729/92 x1^41 x2^42 x3^43"
    ]
    params = ["x1", "x2", "x3"]

    rfl1 = NewRatFunList(3, ratfunlist1)
    rfl2 = ParseRatFun(params, ratfunlist2)
    x = [3057312585776011302, 3795153781312484964, 3415194000889226426]
    if EvaluateRatFunList(rfl1, x, 0) != EvaluateRatFunList(rfl2, x, 0):
        print("- Test failed: functions are different")
        exit(1)
    print("- Test passed!")


def testLSolver(type):
    print("Test {} linear solver".format(type))

    # Solving:
    # [ (a1 + a2) x + (a1 - a2) y == a1/a2,
    #   1/a1 x + 1/a2 y == 0 ]
    # in x,y
    # (the numerical version is the same but sets a1=1, a2=2)
    #
    # This is done by passing this data
    #
    #    (a1 + a2) , (a1 - a2) , a1/a2
    #      1/a1    ,    1/a2   ,  0
    #
    # as a sparse matrix.
    if type == "analytic" or type == "node":
        coeffs = [
            "a1 + a2", "a1 - a2", "a1/a2",
            "1/a1", "1/a2"
        ]
    else:
        coeffs = [
            "3", "-1", "1/2",
            "1", "1/2"
        ]
    cols = [[0,1,2], [0,1]] # non-zero columns for each row
    n_unknowns = 2

    if type == "analytic":
        mygraph, myinput = NewGraphWithInput(len(["a1", "a2"]))
        ccs = ParseIdxRatFun(["a1", "a2"], coeffs)
        sys = AlgAnalyticSparseLSolve(mygraph, myinput, n_unknowns, cols, ccs)
    elif type == "node":
        mygraph, myinput = NewGraphWithInput(len(["a1", "a2"]))
        ccs = ParseRatFun(["a1", "a2"], coeffs)
        rf = AlgRatFunEval(mygraph, myinput, ccs)
        sys = AlgNodeSparseLSolve(mygraph, rf, n_unknowns, cols)
    else:
        mygraph = NewGraph()
        sys = AlgNumericSparseLSolve(mygraph, n_unknowns, cols, coeffs)
    SetOutputNode(mygraph,sys)
    Learn(mygraph)
    if LSolveDepVars(mygraph,sys) != [0,1]:
        print("- Dep vars = ",LSolveDepVars(mygraph,sys))
        print("- Test failed: something wrong with the system")
        exit(1)

    LSolveMarkAndSweepEqs(mygraph, sys)
    if type == "analytic" or type == "numeric":
        LSolveDeleteUnneededEqs(mygraph, sys)

    if type == "analytic" or type == "node":
        rec = ReconstructFunction(mygraph)
        if EvaluateRatFunList(rec, [1,2], 0) != [2767011611056432735, 3689348814741910313]:
            print("- Test failed: something wrong with the reconstructed solution")
            exit(1)
    else:
        if EvaluateGraph(mygraph, [], 0) != [2767011611056432735, 3689348814741910313]:
            print("- Test failed: something wrong with the reconstructed solution")
            exit(1)

    DeleteGraph(mygraph)
    print("- Test passed!")


def testLaurent():
    print("Test Laurent expansion")

    funcs = '(1+x1 x2^2)/(x1+x1 x3^2), (1+x1^2)/(x1+x1^2 x3^2), (x1-x1^3)/(1-x1 x2^2), (1-x1^2)/(x1-x3^2)'.split(',')
    funcs = ParseRatFun(['x1', 'x2', 'x3'], funcs)

    g1,in1 = NewGraphWithInput(3)
    f1 = AlgRatFunEval(g1, in1, funcs)
    SetOutputNode(g1, f1)

    g2,in2 = NewGraphWithInput(2)
    laur = AlgLaurent(g2, in2, g1, 3)
    if LaurentMaxOrders(g2,laur) != [3,3,3,3]:
        print("- Wrong Laurent init")
        exit(1)

    SetOutputNode(g2, laur)
    Learn(g2)
    if LaurentMinOrders(g2,laur) != [-1, -1, 1, 0]:
        print("- Wrong Laurent learn")
        exit(1)

    rec = ReconstructFunction(g2)
    check = [1552403684386711932, 8740806668687644896, 0, 0, 0, 1,
             9223372036854773266, 4100626, 9223372028551007641, 16815129491250,
             1, 576, 331775, 7451573655454030803, 5319074933661716936,
             4183888693711775523, 1696433721680833671]
    if EvaluateRatFunList(rec, [24,45], 10) != check:
        print("- Failed check")
        exit(1)
    DeleteGraph(g1)
    DeleteGraph(g2)
    print("- Test passed")


def testLists():
    print("Test lists")

    list_a = [str(x) for x in range(4)]
    list_b = [str(x) for x in range(5)]
    list_c = [str(x) for x in range(2)]
    list_d = [str(x) for x in range(7)]

    g = NewGraph()
    node_a = AlgRatNumEval(g, list_a)
    node_b = AlgRatNumEval(g, list_b)
    node_c = AlgRatNumEval(g, list_c)
    node_d = AlgRatNumEval(g, list_d)

    node_ab = AlgChain(g, [node_a, node_b])
    node_cd = AlgTake(g, [node_c, node_d],
                      [(0,0), (0,1),
                       (1,0), (1,1), (1,2), (1,3), (1,4), (1,5), (1,6)])

    node_add = AlgAdd(g, [node_ab, node_cd])
    node_mul = AlgMul(g, [node_ab, node_cd])

    node_amat = AlgSlice(g, node_add, 1, 6)
    node_bmat = AlgSlice(g, node_mul, 6)

    mat_mul = AlgSparseMatMul(g, node_amat, node_bmat,
                              2,3,1,
                              [[1,2],[0,1,2]],
                              [[0],[0],[0]])

    tadd = AlgTakeAndAdd(g, [node_amat, node_bmat, mat_mul],
                         [
                             [(0,0),(0,4),(1,2)],
                             [(0,1),(0,2),(2,1)],
                             [(0,3),(2,0)]
                         ])
    SetOutputNode(g, tadd)
    if EvaluateGraph(g, [], 0) != [30, 164, 80]:
        print("- Check failed!")
        exit(1)

    # Test TakeUnique
    unique_ab,fromunique_ab = TakeUnique(g, node_ab)
    SetOutputNode(g, unique_ab)
    if EvaluateGraph(g, [], 0) != list(range(5)) or \
       fromunique_ab != list(range(4)) + list(range(5)):
        print("- Check failed for TakeUnique!")
        exit(1)

    taddbl = AlgTakeAndAddBL(g, [node_amat, node_bmat, mat_mul],
                             [
                                 [(0,0,0,0),(0,4,1,2)],
                                 [(0,1,0,2),(2,1,0,0),(0,1,0,1)],
                                 [(0,3,2,0)]
                            ])
    SetOutputNode(g, taddbl)
    if EvaluateGraph(g, [], 0) != [2*2 + 4*24, 2*4 + 158*2 + 2*2, 2*78]:
        print("- Check failed!")
        exit(1)

    DeleteGraph(g)
    print("- Test passed!")


def testEvaluate():
    print("Testing EvaluatePoints")

    from random import randint
    rand = lambda : randint(123456789123456789, PrimeNo(0)-1)

    funcs = testTutorial2AnalyticCheck()
    g,inp = NewGraphWithInput(3)
    node = AlgRatFunEval(g,inp,ParseRatFun(("x1","x2","x3"),funcs))
    SetOutputNode(g,node)

    points = list(tuple(rand() for _ in range(3)) for _ in range(7))
    EvaluatePoints(g,points)

    DeleteGraph(g)
    print("- Test passed!")


def testUnivariate_(parallel,mod):
    nfuns = 100
    fun = list("1+x^{}".format(i) for i in range(1,nfuns+1))

    g,inp = NewGraphWithInput(1)
    rf = AlgRatFunEval(g,inp,ParseRatFun(["x"],fun))
    SetOutputNode(g,rf)

    if (parallel):
        if (mod):
            res = ParallelReconstructUnivariateMod(g)
        else:
            res = ParallelReconstructUnivariate(g)
    else:
        if (mod):
            res = ReconstructFunctionMod(g)
        else:
            res = ReconstructFunction(g)

    if (mod):
        prime_no = 0
    else:
        prime_no = 10
    pp = PrimeNo(prime_no)
    evals = EvaluateRatFunList(res,[2],prime_no)
    check = list((1+2**i) % pp for i in range(1,nfuns+1))

    if evals == check:
        return SUCCESS
    else:
        print("- Check failed with parallel = {}, mod = {}"
              .format(parallel,mod))
        return ERROR

    DeleteGraph(g)

def testUnivariate():
    print("Testing univariate reconstruction")
    if testUnivariate_(False,False) != SUCCESS:
        exit(1)
    if testUnivariate_(False,True) != SUCCESS:
        exit(1)
    if testUnivariate_(True,False) != SUCCESS:
        exit(1)
    if testUnivariate_(True,True) != SUCCESS:
        exit(1)
    print("- Test passed")

def testNumeric():
    print("Numerical tests")

    # Reconstruction
    length = 10
    ratnums = list('{}/{}'.format(i,i+1) for i in range(1,length))
    g = NewGraph()
    nev = AlgRatNumEval(g, ratnums)
    SetOutputNode(g,nev)
    rec = ReconstructNumeric(g)
    if rec != ratnums:
        print("- Numerical reconstruction failed!")
        exit(1)
    DeleteGraph(g)

    # Chinese remainder
    x = "100000000000000/99999999999"
    y = "100000000000003/99999999998"
    x0 = 8084042683842557234
    x1 = 2691542216189571617
    y0 = 5820300346345761690
    y1 = 4643407741435786205
    cr,ptot = ChineseRemainder([x0,y0],PrimeNo(0),[x1,y1],PrimeNo(1))

    # Rational reconstruction
    rec = RatRec(cr,ptot)
    if rec != [x,y]:
        print("- Test of Chinese remainder + rational reconstruction failed!")
        exit(1)

    # RatMod utility
    if RatMod([x,y], 1) != [x1,y1]:
        print("- Test of RatMod failed!")
        exit(1)

    print("- Test passed")


def testLSolverEx():
    print("Test linear solver 'ex'")
    with GraphContextWithInput(1) as (g,inp):
        unused = AlgRatNumEval(g,["1","2","3"])
        # w0 = -1
        # w1 = 1
        # w2 = t
        # w3 = t^2
        ws = AlgRatFunEval(g,inp,ParseRatFun(["t"],"-1,1,t,t^2".split(",")))
        w1 = (0,1)
        w2 = (0,2)
        w3 = (0,3)
        nonzero_cols = [
            # (w1 * t * x + (w2 * 1/t + w3 * t) * y == (w1))
            [(0,[w1]), (1,[w2,w3]), (2,[w1])],
            # (w1*t - (1)*w2) * y == 0
            [(1,[w1,w2])],
            # (w1 * t * x + (w1*t - (1)*w2) * y == w2)
            [(0,[w1]), (1,[w1,w2]), (2,[w2])]
        ]
        nonzero_ccs = ParseIdxRatFun(
            ["t"],
            ("t,1/t,t,1,"+\
             "t,-1,"+\
             "t,t,-1,1").split(",")
        )
        ls = AlgAnalyticSparseLSolveEx(g,[inp,ws],2,nonzero_cols,nonzero_ccs)
        SetOutputNode(g,ls)
        Learn(g)
        depvars = LSolveDepVars(g, ls)
        indepvars = LSolveIndepVars(g, ls)
        impossible = LSolveIsImpossible(g,ls)
        if impossible or (depvars != [0,1] and depvars != [1,0]):
            print("- Test of depvars failed")
            exit(1)
        LSolveMarkAndSweepEqs(g, ls)
        LSolveDeleteUnneededEqs(g, ls)
        pt = [1234567890]
        rec = ReconstructFunction(g)
        check = ParseRatFun(["t"],["1","(-t+1)/(1+t^3)"])
        if EvaluateRatFunList(rec,pt,0) == EvaluateRatFunList(check,pt,0):
            print("- Test passed")
        else:
            print("- Test failed")
            exit(1)


def testSubgraph():
    print("Test Subgraph Map")
    with GraphContextWithInput(2) as (g,inp):
        funs = ParseRatFun(["x","y"],
                           ["x","x^2","y","y^2"])
        rf = AlgRatFunEval(g,inp,funs)
        SetOutputNode(g,rf)

        with GraphContextWithInput(1) as (g2,inp2):
            rfa = AlgRatFunEval(g2,inp2,ParseRatFun(["x"],["x","x^2"]))
            rfb = AlgRatFunEval(g2,inp2,ParseRatFun(["x"],["x^3","x^4"]))
            sub = AlgSubgraphMap(g2,[rfa,rfb],g)
            SetOutputNode(g2,sub)

            ev = EvaluateGraph(g2,[2],0)
            if ev == [2,2**2,2**2,2**4,2**3,2**6,2**4,2**8]:
                print("- Test passed")
            else:
                print("- Test failed")
                exit(1)


def testSubgraphRec():
    print("Test Subgraph Rec")
    with GraphContextWithInput(3) as (g,inp):
        funs = ParseRatFun(["x","y","t"],
                           ["x + y*t","(1+x*y*t)/(x + y^2*t)"])
        rf = AlgRatFunEval(g,inp,funs)
        SetOutputNode(g,rf)

        with GraphContextWithInput(1) as (g2,inp2):
            sub = AlgSubgraphRec(g2,inp2,g,2)
            SetOutputNode(g2,sub)
            Learn(g2)

            exps = SubgraphRecExponents(g2,sub)
            exps_check = [([(1, 0), (0, 1)], [(0, 0)]),
                          ([(1, 1), (0, 0)], [(0, 2), (1, 0)])]
            if not exps == exps_check:
                print("- Test failed: wrong monomial exponents")

            tval = 1234567890
            ev = EvaluateGraph(g2,[tval],0)
            if ev == [1,tval,1,
                      tval,1,tval,1]:
                print("- Test passed")
            else:
                print("- Test failed")
                exit(1)


def testRatFunEvalFromCoeffs():
    print("Test RatFunEvalFromCoeffs")
    with GraphContextWithInput(2) as (g,inp):

        rn = AlgRatNumEval(g,(1,2,3,4))

        # note: here the coefficients 0,1,3,4 will be indexes inside
        # the first input node (rn) and not actual numerical values.
        fun = ParseRatFun(["x","y"],["(0 x + 1 y)/(2 x + 3 y)"])
        if not len(fun.monomials()[0][0]) == len(fun.monomials()[0][1]) == 2:
            print("- Something wrong with parser function")
            exit(1)

        rf = AlgRatFunEvalFromCoeffs(g,rn,inp,fun)
        SetOutputNode(g,rf)

        rec = ReconstructFunction(g)
        pt = [1234567,67890123]
        check = ParseRatFun(["x","y"],["(1 x + 2 y)/(3 x + 4 y)"])
        if not EvaluateRatFunList(rec,pt,3) == EvaluateRatFunList(check,pt,3):
            print("- Test failed!")
            exit(1)

        print("- Test passed!")


def testRatExpr():
    print("Test RatExprEval")
    rf = RatExprToRatFunList(["(-1/3433683820292512484657849089281 - x^(-1) - y^2)^(-2)",
                              "(x^2 - y^(-1))*(1 + z^(-2))^2"],
                             variables=["x","y","z"])
    point = [758187025378063064,4448666498048535379,2564206275484034988]
    check = [5118394223678577184, 1684865360044594812]
    if EvaluateRatFunList(rf,point,prime_no=0)==check:
        print("- Test passed!")
    else:
        print("- Test failed!")
        exit(1)


if __name__ == '__main__':
    testRatFun()
    testParsing()
    testTutorial2()
    testBasicRatFunInterface()
    testLSolver("analytic")
    testLSolver("node")
    testLSolver("numeric")
    testLSolverEx()
    testLaurent()
    testLists()
    testEvaluate()
    testUnivariate()
    testNumeric()
    testSubgraph()
    testSubgraphRec()
    testRatFunEvalFromCoeffs()
    testRatExpr()
