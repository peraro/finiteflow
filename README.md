FiniteFlow
==========

Installation
------------

### Dependencies

#### CMake

FiniteFlow uses the CMake build system, which comes preinstalled on
many systems, or can be installed with available package managers for
most Linux distrubutions and using Homebrew for macOS.  Binary
installers for several operating systems, as well as additional
information, can be found at [https://cmake.org/](https://cmake.org/).


#### GMP

GMP is preinstalled on many systems.  A `gmp` package is available on
many Linux distrubutions and on Homebrew for macOS.  For more
information see [https://gmplib.org/](https://gmplib.org/).


#### FLINT

FiniteFlow uses a few functions and macros from the FLINT library.
For the installation of FLINT you have two options:

- Install the full FLINT library.  Installation packages are available
  for some Linux distrubutions (including Debian/Ubuntu) and macOS.
  See [http://www.flintlib.org/](http://www.flintlib.org/) for more
  information about FLINT and its installation.

- Alternatively, you can install the `flint-finiteflow-dep` library,
  which is a stripped down version of FLINT which includes the
  dependencies needed by FiniteFlow.  Unlike the full FLINT library,
  this version has no external dependencies besides GMP and can
  therefore be easier to install from source (see the README file of
  `flint-finiteflow-dep` for more information).


### Installing FiniteFlow

FiniteFlow uses the CMake build system.  In order to install
FiniteFlow in a default installation path, you can use the command
```
cmake . && make install
```

In order to install FiniteFlow in a custom path, or specify a
different path for its dependencies use
```
cmake -DCMAKE_INSTALL_PREFIX=/installation/path/prefix \
      -DCMAKE_PREFIX_PATH=/dependencies/installation/path/prefix \
      -DMATHLIBINSTALL=/mathematica/interface/installation/path .
make install
```
where you can omit any of the options in order to pick a default value
for them.

In order to use FiniteFlow from the Mathematica interface, consider
adding the following to your Mathematica `init.m` file
```
(* the next two lines are only needed if a custom MATHLIBINSTALL path was chosen when running cmake *)
$FiniteFlowLibPath = "/mathematica/interface/installation/path"
If[Not[MemberQ[$LibraryPath,$FiniteFlowLibPath]],$LibraryPath = Flatten[{$LibraryPath, $FiniteFlowLibPath }]];

(* the next two lines are needed to locate FiniteFlow.m, if not already in your $Path *)
$FiniteFlowPath = "/path/to/finiteflow/mathlink"
If[Not[MemberQ[$Path,$FiniteFlowPath]],$Path = Flatten[{$Path, $FiniteFlowPath }]];
```
