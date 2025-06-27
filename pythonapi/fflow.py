#! /usr/bin/env python3

'''Module with python bindings for some FiniteFlow routines.

WARNING: This is a w.i.p. and the API is not stable yet.
'''


from _cffi_fflow import lib as _lib, ffi as _ffi
from itertools import chain as _chain
from random import randint as _randint


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

    def _poly_coeffs(self,gettems,getcoeffs,idx):
        nterms = _Check(gettems(self._ptr,idx))
        ccoeffs = getcoeffs(self._ptr,idx)
        ret = list(_ffi.string(ccoeffs[i]).decode() for i in range(nterms))
        _lib.ffFreeCStrArray(ccoeffs)
        return ret

    def num_coeffs(self,idx=None):
        if idx is None:
            return list(self.num_coeffs(i) for i in range(self.size()))
        else:
            return self._poly_coeffs(_lib.ffRatFunNumNTerms,
                                     _lib.ffRatFunNumCoeffs, idx)

    def den_coeffs(self,idx=None):
        if idx is None:
            return list(self.den_coeffs(i) for i in range(self.size()))
        else:
            return self._poly_coeffs(_lib.ffRatFunDenNTerms,
                                     _lib.ffRatFunDenCoeffs, idx)

    def coeffs(self,idx=None):
        if idx is None:
            return list(self.coeffs(i) for i in range(self.size()))
        else:
            return (self.num_coeffs(idx), self.den_coeffs(idx))

    def _poly_exponents(self,getterms,getexps,idx):
        nterms = _Check(getterms(self._ptr,idx))
        nv = self.nvars()
        exps = getexps(self._ptr,idx)
        ret = list(tuple(exps[i*nv: i*nv+nv]) for i in range(nterms))
        _lib.ffFreeMemoryU16(exps)
        return ret

    def num_exponents(self,idx=None):
        if idx is None:
            return list(self.num_exponents(i)
                        for i in range(self.size()))
        else:
            return self._poly_exponents(_lib.ffRatFunNumNTerms,
                                        _lib.ffRatFunNumExponents, idx)

    def den_exponents(self,idx=None):
        if idx is None:
            return list(self.den_exponents(i)
                        for i in range(self.size()))
        else:
            return self._poly_exponents(_lib.ffRatFunDenNTerms,
                                        _lib.ffRatFunDenExponents, idx)

    def exponents(self,idx=None):
        if idx is None:
            return list(self.exponents(i) for i in range(self.size()))
        else:
            return (self.num_exponents(idx), self.den_exponents(idx))

    def num_monomials(self,idx=None):
        if idx is None:
            return list(self.num_monomials(i) for i in range(self.size()))
        else:
            return list(zip(self.num_coeffs(idx), self.num_exponents(idx)))

    def den_monomials(self,idx=None):
        if idx is None:
            return list(self.den_monomials(i) for i in range(self.size()))
        else:
            return list(zip(self.den_coeffs(idx), self.den_exponents(idx)))

    def monomials(self,idx=None):
        if idx is None:
            return list(self.monomials(i) for i in range(self.size()))
        else:
            return (self.num_monomials(idx), self.den_monomials(idx))

    def to_string(self,svars,idx=None):
        if len(svars) != self.nvars():
            raise ValueError("List svars must have the same length as the " +
                             "number of variables self.nvars()")
        if idx is None:
            return list(self.to_string(svars,i) for i in range(self.size()))
        else:
            cvars = [_ffi.new("char[]", x.encode('utf8')) for x in svars]
            cstr = _lib.ffRatFunToStr(self._ptr,idx,cvars)
            if cstr == _ffi.NULL:
                raise ERROR
            ret = _ffi.string(cstr).decode()
            _lib.ffFreeCStr(cstr)
            return ret

class IdxRatFunList:
    '''\
An opaque object representing an indexed list of rational \
functions, namely a list of functions and a list of indexes into \
the list of functions.  The indexes are meant avoid repetition of \
equal entries in the list.

It is used to define sparse linear systems, which often have many \
repeated entries.

As an example using
    functions = [f0,f1,f2]
    indexes = [0,0,1,2,2,0]
would be effectively equivalent to the list
    functions = [f0,f0,f1,f2,f2,f0]
but storing the functions [f0,f1,f2] only once.

Note that the builtin len() function and the .size() method return \
the length of the list of indexes.  Use the method .nfunctions() to \
get the list of unique functions instead.
'''
    def __init__(self):
        self._ptr = _ffi.NULL

    def __del__(self):
        _lib.ffFreeIdxRatFun(self._ptr)

    def size(self):
        if self._ptr == _ffi.NULL:
            return 0
        return _lib.ffIdxRatFunListSize(self._ptr)

    def nfunctions(self):
        if self._ptr == _ffi.NULL:
            return 0
        return _lib.ffIdxRatFunListNFunctions(self._ptr)

    def __len__(self):
        return self.size()

    def nvars(self):
        if self._ptr == _ffi.NULL:
            return 0
        return _lib.ffIdxRatFunListNVars(self._ptr)



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

FFLOW_VERSION = _lib.ffVersion()
FFLOW_VERSION_MINOR = _lib.ffVersionMinor()
__version__ = str(FFLOW_VERSION) + "." + str(FFLOW_VERSION_MINOR)

def MulInv(z, prime_no):
    return _ToUint(_lib.ffMulInv(z, prime_no))

def PrimeNo(i):
    return _ToUint(_lib.ffPrimeNo(i))

def NAvailablePrimes():
    return _lib.ffNAvailablePrimes()

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

def GetOutputNode(graph):
    return _Check(_lib.ffGetOutputNode(graph))

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

def AlgTakeAndAddBL(graph, in_nodes, elems):
    def nested(chain_list):
        for x in chain_list:
            yield x[0]
            yield x[1]
            yield x[2]
            yield x[3]
    return  _Check(_lib.ffAlgTakeAndAddBL(graph, in_nodes, len(in_nodes),
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

def AlgEvalCount(graph, in_node):
    return _Check(_lib.ffAlgEvalCount(graph, in_node))

def EvalCountGet(graph, node):
    return _ToUint(_lib.ffEvalCountGet(graph, node))

def EvalCountReset(graph, node, count=0):
    return _ToUint(_lib.ffEvalCountReset(graph, node, count))


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


def LSolveOnlyNonHomogeneous(graph, node):
    return _StatusCheck(_lib.ffLSolveOnlyNonHomogeneous(graph, node))


def LSolveSparseOutput(graph, node, sparse : True):
    return _StatusCheck(_lib.ffLSolveSparseOutput(graph, node, sparse))


def LSolveSparseOutputWithMaxCol(graph, node, maxcol,
                                 backsubst=True, keepfullout=False):
    return _StatusCheck(_lib.ffLSolveSparseOutputWithMaxCol(graph, node,
                                                            maxcol, backsubst,
                                                            keepfullout))


def LSolveMarkAndSweepEqs(graph, node):
    return _StatusCheck(_lib.ffLSolveMarkAndSweepEqs(graph, node))


def LSolveDeleteUnneededEqs(graph, node):
    return _StatusCheck(_lib.ffLSolveDeleteUnneededEqs(graph, node))


def LSolveIsImpossible(graph,node):
    ret = _Check(_lib.ffLSolveIsImpossible(graph, node))
    return ret == 1


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
    return _StatusCheck(_lib.ffRatFunToJSON(rf._ptr,
                                            json_file.encode('utf8')))


def ParseRatFun(variables, functions):
    '''This uses a simple and limited parser of rational functions:
- functions must be collected under common denominator
- numerators and denominators must be in expanded form
- rational coefficients must be in front of their monomials,
  except for their denominator (e.g. "3/2 z" and "3z/2" are both
  ok but "z 3/2" is not)
This may also return valid rational functions for some invalid
inputs.

If these limitations are too restrictive, consider using the parser of
a proper CAS and then pass the functions to fflow using NewRatFunList
instead.
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


def _getUniqueAndIdx(somelist):
    uni = dict()
    ret = list()
    for el in somelist:
        if el in uni:
            ret.append(uni[el])
        else:
            idx = len(uni)
            ret.append(idx)
            uni[el] = idx
    return (list(k for k,v in sorted(uni.items(), key = lambda x : x[1])), ret)


def ParseIdxRatFun(variables, functions, indexes = None):
    if indexes is None:
        (unique, indexes) = _getUniqueAndIdx(functions)
    else:
        unique = functions
    ret0 = ParseRatFun(variables, unique)
    retc = _lib.ffMoveRatFunToIdx(ret0._ptr, indexes, len(indexes))
    if retc == _ffi.NULL:
        raise ERROR
    ret = IdxRatFunList()
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
        neededlen = len(needed)
    if sum(len(x) for x in non_zero_els) != len(non_zero_coeffs):
        raise ERROR
    retc = _lib.ffAlgAnalyticSparseLSolveIdx(graph, in_node,
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
        neededlen = len(needed)
    if sum(len(x) for x in non_zero_els) != len(non_zero_coeffs):
        raise ERROR
    (uniqueccs, indexes) = _getUniqueAndIdx(non_zero_coeffs)
    coeffs = [_ffi.new("char[]", x.encode('utf8')) for x in uniqueccs]
    retc = _lib.ffAlgNumericSparseLSolve(graph,
                                         len(non_zero_els), n_vars,
                                         [len(x) for x in non_zero_els],
                                         [x for x in _chain(*non_zero_els)],
                                         indexes, coeffs, len(coeffs),
                                         needed, neededlen)
    return _Check(retc)


def AlgNodeSparseLSolve(graph, in_node, n_vars,
                        non_zero_els, needed_vars = None):
    # TODO: add some more checks
    needed = needed_vars
    if needed is None:
        needed = _ffi.NULL
        neededlen = 0
    else:
        neededlen = len(needed)
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
                            + " for both the numerator and the denominator.")
        n_num_terms.append(len(ratfun[0]))
        n_den_terms.append(len(ratfun[1]))
        for poly in ratfun:
            for t in poly:
                coeffs.append(_ffi.new("char[]", t[0].encode('utf8')))
                if len(t[1]) != nvars:
                    raise TypeError("The length of the exponents' list is not "
                                    + str(nvars))
                exponents.extend(t[1])
    resc = _lib.ffNewRatFunList(nvars, len(allterms), n_num_terms, n_den_terms,
                                coeffs, exponents)
    if resc == _ffi.NULL:
        raise ERROR
    res = RatFunList()
    res._ptr = resc
    return res


def NewIdxRatFunList(nvars, allterms, indexes = None):
    if indexes is None:
        (unique, indexes) = _getUniqueAndIdx(allterms)
    else:
        unique = allterms
    ret0 = NewRatFunList(nvars, unique)
    retc = _lib.ffMoveRatFunToIdx(ret0._ptr, indexes, len(indexes))
    if retc == _ffi.NULL:
        raise ERROR
    ret = IdxRatFunList()
    ret._ptr = retc
    return ret


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
    '''\
ReconstructFunction(graph,**kwargs) reconstructs the rational function
defined by the graph. The allowed keyword arguments kwargs are:
start_mod, min_primes, max_primes, max_deg, dbginfo, polymethod,
n_threads.  All of them must have integer values.
'''
    recopt = _ffi.new("FFRecOptions *",kwargs)
    res = _ffi.new("FFRatFunList **")
    ret = _to_status[_lib.ffReconstructFunction(graph,recopt[0],res)]
    if ret is SUCCESS:
        ret = RatFunList()
        ret._ptr = res[0]
    return ret

def ReconstructFunctionMod(graph, **kwargs):
    '''\
ReconstructFunctionMod(graph,**kwargs) reconstructs the rational
function defined by the graph modulo the prime PrimeNo(start_mod),
with start_mod=0 by default. The allowed keyword arguments kwargs are:
start_mod, max_deg, dbginfo, polymethod, n_threads.  All of them must
have integer values.
'''
    recopt = _ffi.new("FFRecOptions *",kwargs)
    res = _ffi.new("FFRatFunList **")
    ret = _to_status[_lib.ffReconstructFunctionMod(graph,recopt[0],res)]
    if ret is SUCCESS:
        ret = RatFunList()
        ret._ptr = res[0]
    return ret

def ReconstructFromCurrentEvaluations(graph, **kwargs):
    recopt = _ffi.new("FFRecOptions *",kwargs)
    res = _ffi.new("FFRatFunList **")
    ret = _to_status[_lib.ffReconstructFromCurrentEvaluations(graph,recopt[0],res)]
    if ret is SUCCESS:
        ret = RatFunList()
        ret._ptr = res[0]
    return ret

def ReconstructFromCurrentEvaluationsMod(graph, **kwargs):
    recopt = _ffi.new("FFRecOptions *",kwargs)
    res = _ffi.new("FFRatFunList **")
    ret = _to_status[_lib.ffReconstructFromCurrentEvaluationsMod(graph,recopt[0],res)]
    if ret is SUCCESS:
        ret = RatFunList()
        ret._ptr = res[0]
    return ret

def AllDegrees(graph,**kwargs):
    recopt = _ffi.new("FFRecOptions *",kwargs)
    degptr = _lib.ffAllDegrees(graph,recopt[0])
    if degptr == _ffi.NULL:
        raise ERROR
    nout = GraphNParsOut(graph)
    ret = iter(_ffi.unpack(degptr, 2*nout))
    _lib.ffFreeMemoryU32(degptr)
    return [(next(ret),next(ret)) for i in range(nout)]

def DumpDegrees(graph,filename):
    cfile = _ffi.new("char[]", filename.encode('utf8'))
    return _StatusCheck(_lib.ffDumpDegrees(graph,cfile))

def NParsFromDegreeFile(filename):
    nparsin = _ffi.new('unsigned[1]')
    nparsout = _ffi.new('unsigned[1]')
    cfile = _ffi.new("char[]", filename.encode('utf8'))
    ret = _StatusCheck(_lib.ffNParsFromDegreeFile(cfile,nparsin,nparsout))
    return (nparsin[0],nparsout[0])

def LoadDegrees(graph,filename):
    cfile = _ffi.new("char[]", filename.encode('utf8'))
    return _StatusCheck(_lib.ffLoadDegrees(graph,cfile))

def LoadEvaluations(graph,files):
    cfiles = [_ffi.new("char[]", f.encode('utf8')) for f in files]
    return _StatusCheck(_lib.ffLoadEvaluations(graph,cfiles,len(files)))

def DumpSamplePoints(graph,filename,**kwargs):
    recopt = _ffi.new("FFRecOptions *",kwargs)
    cfile = _ffi.new("char[]", filename.encode('utf8'))
    return _StatusCheck(_lib.ffDumpSamplePoints(graph,cfile,recopt[0]))

def NSamplePointsInFile(filename):
    cfile = _ffi.new("char[]", filename.encode('utf8'))
    return _Check(_lib.ffNSamplePointsInFile(cfile))

def EvaluatePointsInFile(graph,filename,start,n_points,n_threads=0):
    cfile = _ffi.new("char[]", filename.encode('utf8'))
    return _StatusCheck(_lib.ffEvaluatePointsInFile(graph,cfile,
                                                    start,n_points,n_threads))

def DumpEvaluations(graph,filename):
    cfile = _ffi.new("char[]", filename.encode('utf8'))
    return _StatusCheck(_lib.ffDumpEvaluations(graph,cfile))

def EvaluatePoints(graph,points,prime_no=0,n_threads=0):

    nparsin = NodeNParsOut(graph,0)

    def append_tuple(x, el):
        return x + (el,)
    def append_list(x, el):
        return x + [el]
    app_ = {list : append_list, tuple : append_tuple}
    appto = lambda x : x if len(x) == nparsin+1 else app_[type(x)](x,prime_no)

    full_pts = list(_chain(*map(appto,points)))
    ret = _lib.ffEvaluatePoints(graph,full_pts,len(points),n_threads)
    if ret == _ffi.NULL:
        raise ERROR

    nparsout = GraphNParsOut(graph)
    res = list(tuple(ret[i:i+nparsout]) for i in range(len(points)))
    _lib.ffFreeMemoryU64(ret)

    return res


def TakeUnique(graph, nodein, nevals=3):
    '''\
TakeUnique(graph, nodein, nevals=3) returns a tuple

    (newnode, fromouts)

where `newnode` is the integer id of a newly created node which takes
`nodein` as input and returns the subset of its unique elements.
These are identified via numerical evaluations (their number can be
changed with the optional argument `nevals`).  The list `fromouts` is
such that
```
    input == list(output[i] for i in fromouts)
```
where input/output is the input/output of the new node.
'''
    SetOutputNode(graph, nodein)
    nout = GraphNParsOut(graph)
    nin = NodeNParsOut(graph,0)
    rand = lambda : _randint(123456789123456789, PrimeNo(200)-1)
    evals = list(EvaluateGraph(graph, list(rand() for _ in range(nin)), 0)
                 for _ in range(nevals))
    evals = list(map(tuple, zip(*evals))) # transpose
    asso = dict()
    outn = 0
    outs = []
    fromouts = []
    for ii in range(nout):
        evii = evals[ii]
        if not evii in asso:
            asso[evii] = outn
            outs.append(ii)
            fromouts.append(outn)
            outn += 1
        else:
            fromouts.append(asso[evii])
    newnode = AlgTake(graph, [nodein], list((0,jj) for jj in outs))
    return (newnode, fromouts)


if __name__ == '__main__':
    pass
