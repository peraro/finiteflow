#! /usr/bin/env python3

'''Module with python bindings for some FiniteFlow routines.

WARNING: This is a w.i.p. and the API is not stable yet.
'''


from _cffi_fflow import lib as _lib, ffi as _ffi
from itertools import chain as _chain


_FF_SUCCESS = 0
_FF_ERROR = 2**32-1
_FF_MISSING_POINTS = _FF_ERROR - 1
_FF_MISSING_PRIMES = _FF_ERROR - 2
_FF_FAILED = 2**64 - 1


class Singleton(type):
    _instances = {}
    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]


class _SuccessType(metaclass=Singleton):
    def __str__(self):
        return 'SUCCESS'
    def __repr__(self):
        return 'fflow.SUCCESS'


class _ErrorType(Exception,metaclass=Singleton):
    def __str__(self):
        return 'ERROR'
    def __repr__(self):
        return 'fflow.ERROR'


class _MissingPointsType(metaclass=Singleton):
    def __str__(self):
        return 'MISSING_POINTS'
    def __repr__(self):
        return 'fflow.MISSING_POINTS'


class _MissingPrimesType(metaclass=Singleton):
    def __str__(self):
        return 'MISSING_PRIMES'
    def __repr__(self):
        return 'fflow.MISSING_PRIMES'


class _FailedType(Exception,metaclass=Singleton):
    def __str__(self):
        return 'FAILED'
    def __repr__(self):
        return 'fflow.FAILED'


class RatFunList:
    def __init__(self):
        self._ptr = _ffi.NULL

    def __del__(self):
        _lib.ffFreeRatFun(self._ptr)

    def size(self):
        if self._ptr == _ffi.NULL:
            return 0
        return _lib.ffRatFunListSize(self._ptr)

    def __len__(self):
        return self.size()

    def nvars(self):
        if self._ptr == _ffi.NULL:
            return 0
        return _lib.ffRatFunListNVars(self._ptr)



SUCCESS = _SuccessType()
ERROR = _ErrorType()
MISSING_POINTS = _MissingPointsType()
MISSING_PRIMES = _MissingPrimesType()
FAILED = _FailedType()


_to_status = {
    _FF_SUCCESS : SUCCESS,
    _FF_ERROR : ERROR,
    _FF_MISSING_PRIMES : MISSING_PRIMES,
    _FF_MISSING_POINTS : MISSING_POINTS
}

def _ToUint(z):
    if z == _FF_FAILED:
        raise FAILED
    return z

def _Check(arg):
    if _lib.ffIsError(arg):
        raise _to_status[arg]
    return arg

def _CheckIter(arg):
    for v in arg:
        if _lib.ffIsError(v):
            raise to_status[v]
    return arg

def _StatusCheck(arg):
    if _lib.ffIsError(arg):
        raise ERROR
    return _to_status[arg]


_lib.ffInit()

def MulInv(z, prime_no):
    return _ToUint(_lib.ffMulInv(z, prime_no))

def PrimeNo(i):
    return _ToUint(_lib.PrimeNo(i))

def DefaultNThreads():
    return _lib.ffDefaultNThreads()

def NewGraph():
    return _lib.ffNewGraph()

def NewGraphWithInput(nvars):
    innode = _ffi.new('FFStatus[1]')
    graph = _lib.ffNewGraphWithInput(nvars, innode)
    return _CheckIter((graph, innode[0]))

def NewGraphDummy(n_in, n_out):
    return _lib.ffNewGraphDummy(n_in, n_out)

def DeleteGraph(graph):
    return _StatusCheck(_lib.ffDeleteGraph(graph))

def DeleteNode(graph, node):
    return _StatusCheck(_lib.ffDeleteNode(graph,node))

def SetOutputNode(graph, node):
    return _StatusCheck(_lib.ffSetOutputNode(graph,node))

def SetGraphInput(graph, nvars):
    return _Check(_lib.ffSetGraphInput(graph,nvars))

def GraphNParsOut(graph):
    return _Check(_lib.ffGraphNParsOut(graph))

def NodeNParsOut(graph, node):
    return _Check(_lib.ffNodeNParsOut(graph,node))

def MakeNodeMutable(graph, node):
    return _StatusCheck(_lib.ffMakeNodeMutable(graph,node))

def PruneGraph(graph):
    return _StatusCheck(_lib.ffPruneGraph(graph))

def Learn(graph):
    return _StatusCheck(_lib.ffLearn(graph))

def EvaluateGraph(graph,z,primeno):
    nparsin = _lib.ffNodeNParsOut(graph,0)
    if _lib.ffIsError(nparsin) or len(z) != nparsin:
        raise FAILED
    retc = _lib.ffEvaluateGraph(graph,z,primeno)
    if retc == _ffi.NULL:
        raise FAILED
    ret = _ffi.unpack(retc,_lib.ffGraphNParsOut(graph))
    _lib.ffFreeMemoryU64(retc)
    return ret

def AlgSimpleSubgraph(graph, in_nodes, subgraph):
    return _Check(_lib.ffAlgSimpleSubgraph(graph, in_nodes, len(in_nodes),
                                           subgraph))

def AlgMemoizedSubgraph(graph, in_nodes, subgraph):
    return _Check(_lib.ffAlgMemoizedSubgraph(graph,in_nodes,len(in_nodes),
                                             subgraph))

def AlgJSONSparseLSolve(graph, in_node, json_file):
    return _Check(_lib.ffAlgJSONSparseLSolve(graph, in_node,
                                             json_file.encode('utf8')))

def AlgJSONRatFunEval(graph, in_node, json_file):
    return _Check(_lib.ffAlgJSONRatFunEval(graph, in_node,
                                           json_file.encode('utf8')))

def AlgRatFunEval(graph, in_node, rf):
    if not type(rf) is RatFunList:
        raise TypeError("Third argument of AlgRatFunEval() must be a RatFunList")
    return _Check(_lib.ffAlgRatFunEval(graph, in_node, rf._ptr))

def AlgRatNumEval(graph, nums):
    cnums = [_ffi.new("char[]", x.encode('utf8')) for x in nums]
    return _Check(_lib.ffAlgRatNumEval(graph, cnums, len(nums)))

def AlgLaurent(graph, in_node, subgraph, order, max_deg=-1):
    if type(order) is list:
        nout = _lib.ffGraphNParsOut(subgraph)
        if _lib.ffIsError(nout) or len(z) != nout:
            raise ERROR
        return _Check(_lib.ffAlgLaurent(graph, in_node, subgraph,
                                       order, max_deg))
    else:
        return _Check(_lib.ffAlgLaurentConstOrder(graph, in_node, subgraph,
                                                  order, max_deg))

def AlgMatMul(graph, in_node_a, in_node_b, n_rows_a, n_cols_a, n_cols_b):
    return _Check(_lib.ffAlgMatMul(graph, in_node_a, in_node_b,
                                   n_rows_a, n_cols_a, n_cols_b))

def AlgChain(graph, in_nodes):
    return _Check(_lib.ffAlgChain(graph, in_nodes, len(in_nodes)))

def AlgTake(graph, in_nodes, elems):
    for x in elems:
        if len(x) != 2:
            raise ERROR
    return _Check(_lib.ffAlgTake(graph, in_nodes, len(in_nodes),
                                 [x for x in _chain(*elems)],
                                 len(elems)))

def AlgSlice(graph, in_node, begin, end=-1):
    return _Check(_lib.ffAlgSlice(graph, in_node, begin, end))

def AlgAdd(graph, in_nodes):
    return _Check(_lib.ffAlgAdd(graph, in_nodes, len(in_nodes)))

def AlgMul(graph, in_nodes):
    return _Check(_lib.ffAlgMul(graph, in_nodes, len(in_nodes)))

def AlgTakeAndAdd(graph, in_nodes, elems):
    def nested(chain_list):
        for x in chain_list:
            yield x[0]
            yield x[1]
    return  _Check(_lib.ffAlgTakeAndAdd(graph, in_nodes, len(in_nodes),
                                        len(elems),
                                        [len(x) for x in elems],
                                        [y for y in nested(_chain(*elems))]))

def AlgSparseMatMul(graph, in_node_a, in_node_b, n_rows_a, n_cols_a, n_cols_b,
                    non_zeroes_a, non_zeroes_b):
    if len(non_zeroes_a) != n_rows_a or len(non_zeroes_b) != n_cols_a:
        raise ERROR
    return _Check(_lib.ffAlgSparseMatMul(graph, in_node_a, in_node_b,
                                         n_rows_a, n_cols_a, n_cols_b,
                                         [len(x) for x in non_zeroes_a],
                                         [x for x in _chain(*non_zeroes_a)],
                                         [len(x) for x in non_zeroes_b],
                                         [x for x in _chain(*non_zeroes_b)]))


def LaurentMaxOrders(graph, node):
    retc = _lib.ffLaurentMaxOrders(graph, node)
    if retc == _ffi.NULL:
        raise ERROR
    res =  _ffi.unpack(retc,_lib.ffSubgraphNParsout(graph, node))
    _lib.ffFreeMemoryS32(retc)
    return res

def LaurentMinOrders(graph, node):
    retc = _lib.ffLaurentMinOrders(graph, node)
    if retc == _ffi.NULL:
        raise ERROR
    res =  _ffi.unpack(retc,_lib.ffSubgraphNParsout(graph, node))
    _lib.ffFreeMemoryS32(retc)
    return res


def LSolveResetNeededVars(graph, node, needed_vars):
    return _StatusCheck(_lib.ffLSolveResetNeededVars(graph, node, needed_vars, len(needed_vars)))


def LSolveOnlyHomogeneous(graph, node):
    return _StatusCheck(_lib.ffLSolveOnlyHomogeneous(graph, node))


def LSolveSparseOutput(graph, node, sparse : True):
    return _StatusCheck(_lib.ffLSolveSparseOutput(graph, node, sparse))


def LSolveMarkAndSweepEqs(graph, node):
    return _StatusCheck(_lib.ffLSolveMarkAndSweepEqs(graph, node))


def LSolveDeleteUnneededEqs(graph, node):
    return _StatusCheck(_lib.ffLSolveDeleteUnneededEqs(graph, node))


def LSolveNDepVars(graph, node):
    return _Check(_lib.ffLSolveNDepVars(graph, node))


def LSolveDepVars(graph, node):
    retc = _lib.ffLSolveDepVars(graph, node)
    if retc == _ffi.NULL:
        raise ERROR
    ret = _ffi.unpack(retc, _lib.ffLSolveNDepVars(graph, node))
    _lib.ffFreeMemoryU32(retc)
    return ret


def LSolveNIndepVars(graph, node, i=0):
    return _Check(_lib.ffLSolveNIndepVars(graph, node, i))


def LSolveIndepVars(graph, node, i=0):
    retc = _lib.ffLSolveIndepVars(graph, node, i)
    if retc == _ffi.NULL:
        raise ERROR
    ret = _ffi.unpack(retc, _lib.ffLSolveNIndepVars(graph, node, i))
    _lib.ffFreeMemoryU32(retc)
    return ret


def RatFunToJSON(rf, json_file):
    if not type(rf) is RatFunList:
        raise TypeError("First argument of RatFunToJSON() must be a RatFunList")
    return _StatusCheck(_lib.ffEvaluateRatFunList(rf._ptr,
                                                  json_file.encode('utf8')))


def ParseRatFun(variables, functions):
    '''This uses a simple and limited parser of rational functions:
- functions must be collected under common denominator
- numerators and denominators must be in expanded form
- rational coefficients must be in front of their monomials
  (e.g. "1/2 z" is ok but "z/2" is not)

If these lmitations are too restrictive, consider using the parser of
a proper CAS and then pass the functions to fflow using
ffNewRatFunList instead.
    '''
    rf = [_ffi.new("char[]", x.encode('utf8')) for x in functions]
    rfl = [len(x) for x in functions]
    z = [_ffi.new("char[]", x.encode('utf8')) for x in variables]

    retc = _lib.ffParseRatFunEx(z, len(z), rf, rfl, len(rf))
    if retc == _ffi.NULL:
        raise ERROR
    ret = RatFunList()
    ret._ptr = retc
    return ret


def AlgAnalyticSparseLSolve(graph, in_node, n_vars,
                            non_zero_els, non_zero_coeffs,
                            needed_vars = None):
    needed = needed_vars
    if needed is None:
        needed = _ffi.NULL
        neededlen = 0
    else:
        needelen = len(needed)
    if sum(len(x) for x in non_zero_els) != len(non_zero_coeffs):
        raise ERROR
    retc = _lib.ffAlgAnalyticSparseLSolve(graph, in_node,
                                          len(non_zero_els), n_vars,
                                          [len(x) for x in non_zero_els],
                                          [x for x in _chain(*non_zero_els)],
                                          non_zero_coeffs._ptr,
                                          needed,
                                          neededlen)
    return _Check(retc)


def AlgNumericSparseLSolve(graph, n_vars,
                           non_zero_els, non_zero_coeffs,
                           needed_vars = None):
    needed = needed_vars
    if needed is None:
        needed = _ffi.NULL
        neededlen = 0
    else:
        needelen = len(needed)
    if sum(len(x) for x in non_zero_els) != len(non_zero_coeffs):
        raise ERROR
    coeffs = [_ffi.new("char[]", x.encode('utf8')) for x in non_zero_coeffs]
    retc = _lib.ffAlgNumericSparseLSolve(graph,
                                         len(non_zero_els), n_vars,
                                         [len(x) for x in non_zero_els],
                                         [x for x in _chain(*non_zero_els)],
                                         coeffs,
                                         needed,
                                         neededlen)
    return _Check(retc)


def AlgNodeSparseLSolve(graph, in_node, n_vars,
                        non_zero_els, needed_vars = None):
    # TODO: add some more checks
    needed = needed_vars
    if needed is None:
        needed = _ffi.NULL
        neededlen = 0
    else:
        needelen = len(needed)
    retc = _lib.ffAlgNodeSparseLSolve(graph, in_node,
                                      len(non_zero_els), n_vars,
                                      [len(x) for x in non_zero_els],
                                      [x for x in _chain(*non_zero_els)],
                                      needed,
                                      neededlen)
    return _Check(retc)


def NewRatFunList(nvars, allterms):
    coeffs = []
    exponents = []
    n_num_terms = []
    n_den_terms = []
    for ratfun in allterms:
        if len(ratfun) != 2:
            raise TypeError("For each rational function, terms should be listed"
                            + "for both the numerator and the denominator.")
        n_num_terms.append(len(ratfun[0]))
        n_den_terms.append(len(ratfun[1]))
        for poly in ratfun:
            for t in poly:
                coeffs.append(_ffi.new("char[]", t[0].encode('utf8')))
                if len(t[1]) != nvars:
                    raise TypeError("The list of exponents the length is not "
                                    + str(nvars))
                exponents.extend(t[1])
    resc = _lib.ffNewRatFunList(nvars, len(allterms), n_num_terms, n_den_terms,
                                coeffs, exponents)
    if resc == _ffi.NULL:
        raise ERROR
    res = RatFunList()
    res._ptr = resc
    return res


def EvaluateRatFunList(rf, z, prime_no):
    if not type(rf) is RatFunList:
        raise TypeError("First argument of EvaluateRatFunList() must be a RatFunList")
    if rf._ptr == _ffi.NULL or rf.nvars() != len(z):
        raise FAILED
    retc = _lib.ffEvaluateRatFunList(rf._ptr, z, prime_no)
    if retc == _ffi.NULL:
        raise FAILED
    ret = _ffi.unpack(retc, _lib.ffRatFunListSize(rf._ptr))
    _lib.ffFreeMemoryU64(retc)
    return ret


def ReconstructFunction(graph, **kwargs):
    recopt = _ffi.new("FFRecOptions *",kwargs)
    res = _ffi.new("FFRatFunList **")
    ret = _to_status[_lib.ffReconstructFunction(graph,recopt[0],res)]
    if ret is SUCCESS:
        ret = RatFunList()
        ret._ptr = res[0]
    return ret


if __name__ == '__main__':
    pass
