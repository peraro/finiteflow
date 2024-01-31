#! /usr/bin/env python3

import pathlib
from cffi import FFI

thisfile = pathlib.Path(__file__).resolve()

def getFFlowPath():
    try:
        with open(str(thisfile.parent.parent / 'install_manifest.txt')) as f:
            for l in f:
                path = pathlib.Path(l)
                if "libfflow" in path.name:
                    return path.resolve()
    except:
        pass
    raise RuntimeError("FiniteFlow library not found")

ffibuilder = FFI()

fflowpath = getFFlowPath()
header = open(str((thisfile.parent.parent / 'include' / 'fflow' / 'capi.h').resolve())).read()
header = header.split("/* API begin */")[1].split("/* API end */")[0]

ffibuilder.cdef(header)

fflowlibdir = str(fflowpath.parent.resolve())
includedir = str(thisfile.parent.parent / 'include')

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
