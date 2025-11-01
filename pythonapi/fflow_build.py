#! /usr/bin/env python3

import pathlib
from cffi import FFI
import sys

thisfile = pathlib.Path(__file__).resolve()
only_source = len(sys.argv)>=2 and sys.argv[1] == '--source-only'

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

if not only_source:
    fflowpath = getFFlowPath()
header = open(str((thisfile.parent.parent / 'include' / 'fflow' / 'capi.h').resolve())).read()
header = header.split("/* API begin */")[1].split("/* API end */")[0]

ffibuilder.cdef(header)

if not only_source:
    fflowlibdir = str(fflowpath.parent.resolve())
includedir = str(thisfile.parent.parent / 'include')

if only_source:
    extra_link_args = []
    libraries = []
else:
    extra_link_args = ["-L" + fflowlibdir, "-Wl,-rpath," + fflowlibdir]
    libraries = ['fflow']

ffibuilder.set_source("_cffi_fflow",
                      r'''
                      #include <fflow/capi.h>
                      ''',
                      libraries = libraries,
                      include_dirs=[includedir],
                      verbose=True,
                      extra_link_args=extra_link_args
                      )

if __name__ == "__main__":
    pathlib.Path("_cffi_fflow.c").touch()
    if not only_source:
        ffibuilder.compile(verbose=True)
    else:
        ffibuilder.emit_c_code("_cffi_fflow.c")
