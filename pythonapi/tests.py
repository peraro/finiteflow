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
        print("- Test 1/2 failed :(")
        exit(1)
    print("- Test 1/2 passed!!!")

    # build another graph evaluating the reconstructed result
    graph2,innode2 = NewGraphWithInput(nvars)
    rf2 = AlgRatFunEval(graph2, innode2, rec)
    SetOutputNode(graph2, rf2)
    output3 = EvaluateGraph(graph2, x, prime_no)

    # compare again
    if output1 != output3:
        print("- Test 2/2 failed :(")
        exit(1)
    print("- Test 2/2 passed!!!")


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
        ccs = ParseRatFun(["a1", "a2"], coeffs)
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
        print(- "Test failed: something wrong with the system")
        exit(1)

    LSolveMarkAndSweepEqs(mygraph, sys)
    if type == "analytic" or type == "numeric":
        LSolveDeleteUnneededEqs(mygraph, sys)

    if type == "analytic" or type == "node":
        rec = ReconstructFunction(mygraph)
        if EvaluateRatFunList(rec, [1,2], 0) != [2767011611056432735, 3689348814741910313]:
            print(- "Test failed: something wrong with the reconstructed solution")
            exit(1)
    else:
        if EvaluateGraph(mygraph, [], 0) != [2767011611056432735, 3689348814741910313]:
            print(- "Test failed: something wrong with the reconstructed solution")
            exit(1)

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
    print("- Test passed!")



if __name__ == '__main__':
    testRatFun()
    testParsing()
    testTutorial2()
    testBasicRatFunInterface()
    testLSolver("analytic")
    testLSolver("node")
    testLSolver("numeric")
    testLaurent()
    testLists()
