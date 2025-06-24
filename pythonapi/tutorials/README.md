Tutorials using the Python API
==============================

This is a collection of tutorials using the Python API of the
FiniteFlow C/C++ library. This API covers a significant subset of the
Mathematica one.  Note that this API does not rely on a computer
algebra system, hence multi-precision numbers are passed as
strings. Rational functions can be passed as monomial data or as
strings to be parsed, although the parser of FiniteFlow has
limitations.

We provide the following tutorials, which are Python ports of those
included in the [Mathematica one](../../mathlink/tutorial.wl):
- [`1_basic_usage.py`](1_basic_usage.py) is a minimal introduction to
  the package which describes a very basic basic example.
- [`2_linear_systems.py`](2_linear_systems.py) shows how to use an
  algorithm with a *learning phase*, more specifically a linear solver.
- [`3_laurent.py`](3_laurent.py) shows how to use an example using a
  *subgraph* to compute Laurent expansions.

We recommend to read them in order.
