FAQ
===

The installation fails with the error "`cannot find -luuid`"
----------------------------------------------------------

This is generally a MathLink issue (not specific to FiniteFlow) on
some Linux distributions.  It can be solved by installing the
`uuid-dev` package or an equivalent one, depending on the
distribution.


A procedure fails, saying the arguments are not polynomials or rational functions
---------------------------------------------------------------------------------

When passing rational functions to dataflow procedures, they should
be collected under a common denominator, i.e. written as ratios of
two polynomials.  This is also true for the coefficients multiplying
the unknowns of linear systems, where one can use the option
`"ApplyFunction"` to fix the issue.  As an example, this will cause an
error
```
  FFDenseSolve[{(t + 1/t) x + y == 3}, {x,y}]
```
because the coefficient of `x` is not written as a ratio of two
polynomials.  However the following works
```
  FFDenseSolve[{(t + 1/t) x + y == 3}, {x,y}, "ApplyFunction"->Together]
```
and is equivalent to
```
  FFDenseSolve[{Together[t + 1/t] x + y == 3}, {x,y}]
```


The definition of a node in a graph returns `$Failed`
-----------------------------------------------------

This can happen because of multiple reasons, but it is often due to the fact that the expected numbers of inputs of the new node, or their length, is incompatible with the selected input nodes.
