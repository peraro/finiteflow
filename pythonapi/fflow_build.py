#! /usr/bin/env python3

import pathlib
from cffi import FFI
from ctypes.util import find_library

fflowpath = find_library("fflow")
if fflowpath is None:
    raise RuntimeError("FiniteFlow library not found")

ffibuilder = FFI()

header = open(str((pathlib.Path(__file__).parent.parent / 'include' / 'fflow' / 'capi.h').resolve())).read()
header = header.split("/* API begin */")[1].split("/* API end */")[0]

ffibuilder.cdef(header)

fflowlibdir = str(pathlib.Path(fflowpath).parent.resolve())
includedir = str(pathlib.Path(__file__).parent.parent / 'include')

ffibuilder.set_source("_cffi_fflow",
                      r'''
                      #include <fflow/capi.h>
                      ''',
                      libraries = ['fflow'],
                      include_dirs=[includedir],
                      extra_link_args=["-L" + fflowlibdir,
                                       "-Wl,-rpath," + fflowlibdir]
                      )

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
