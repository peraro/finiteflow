FAQ
===

The installation fails with the error "`cannot find -luuid`"
----------------------------------------------------------

This is generally a MathLink issue (not specific to FiniteFlow) on
some Linux distributions.  It can be solved by installing the
`uuid-dev` package or an equivalent one, depending on the
distribution.


The installation procedure is unsuccessful
------------------------------------------
If you encounter issues with the installation, including failure to load the library or the Mathematica package, please follow the steps below before contacting the developer.  These may help you identifying the problem or provide useful diagnostic information to include in support requests.
- Clean up the current installation using `make clean && make uninstall` followed by `rm CMakeCache.txt`.
- Update FiniteFlow to the latest version and repeat the installation procedure, carefully checking all specified paths (if any) are correct.  Monitor the output of `cmake` and `make` for errors.  If errors occur that you are unable to resolve, report them to the developer.
- If you plan to use the Mathematica interface, make sure that Mathematica is found during the `cmake` configuration and that its version matches the one you intend to use.  If not, define the `Mathematica_FIND_VERSION` and/or `Mathematica_FIND_VERSION_EXACT` variables (e.g., with `cmake -DMathematica_FIND_VERSION="13.2.1" -DMathematica_FIND_VERSION_EXACT=TRUE ...`) to specify the required Mathematica version.  Refer to the [FindMathematica manual](https://github.com/sakra/FindMathematica/blob/master/MANUAL.md#used-variables) for more information.
- Execute the command `make tests && ./testjson`.  If the output ends with `Usage: ./testjson input output` then the C++ installation was successful.  Otherwise, check your installation procedure again and, if the errors persist, report the outcome to the developer.
- Execute the command `math -script mathlink/tests.m`.  You may need to replace `math` with the command for executing Mathematica or the Wolfram kernel from command line (e.g. `wolframscript`).  If the tests pass without any error, then the Mathematica interface has been installed correctly.
- If the previous command reported failure in loading the FiniteFlow library, execute `LD_DEBUG=libs math -script mathlink/tests.m` which repeats the previous test but may report more information about this problem.  This information is often helpful in detecting issues, such having multiple conflicting versions of FiniteFlow or its dependencies installed.  If you're still unable to resolve the issue, include the output of this command when reporting the problem to the developer.


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
