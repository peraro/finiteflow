'''Updates the embedded files of the fmt library taking them from
github and replacing their namespace with fflow.

'''

import urllib2, re, os


def is_header(fname):
    "Check if file is a header"
    return fname.endswith(".h")


def main():
    "Main function"

    join_path = os.path.join
    thisdir = os.path.dirname(os.path.realpath(__file__))
    builddir = os.path.abspath(join_path(thisdir, os.pardir))

    regexes = [(r'#include "format\.h"', '#include <fflow/format.h>'),
               (r'#include "ostream\.h"', '#include <fflow/ostream.h>'),
               (r'\bfmt\b', 'fflow')]

    fmtbasicurl = "https://raw.githubusercontent.com/fmtlib/fmt/master/fmt/{}"
    fmtfiles = ["format.h", "format.cc", "ostream.h", "ostream.cc"]

    headerdir = join_path(builddir, 'include', 'fflow')
    srcdir = join_path(builddir, 'src')

    for fname in fmtfiles:
        destdir = headerdir if is_header(fname) else srcdir
        src = urllib2.urlopen(fmtbasicurl.format(fname)).read()
        for regex in regexes:
            src = re.sub(regex[0], regex[1], src)
        with open(join_path(destdir, fname), 'w') as dest:
            dest.write(src)


if __name__ == "__main__":
    main()
