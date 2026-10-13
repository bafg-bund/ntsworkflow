# Set of functions in example RcppArmadillo package

These four functions are created when `RcppArmadillo.package.skeleton()`
is invoked to create a skeleton packages.

## Usage

``` r
rcpparma_hello_world()
rcpparma_outerproduct(x)
rcpparma_innerproduct(x)
rcpparma_bothproducts(x)
```

## Arguments

- x:

  a numeric vector

## Value

`rcpparma_hello_world()` does not return a value, but displays a message
to the console.

`rcpparma_outerproduct()` returns a numeric matrix computed as the outer
(vector) product of `x`.

`rcpparma_innerproduct()` returns a double computer as the inner
(vector) product of `x`.

`rcpparma_bothproducts()` returns a list with both the outer and inner
products.

## Details

These are example functions which should be largely self-explanatory.
Their main benefit is to demonstrate how to write a function using the
Armadillo C++ classes, and to have to such a function accessible from R.

## References

See the documentation for Armadillo, and RcppArmadillo, for more
details.

## Author

Dirk Eddelbuettel

## Examples

``` r
  x <- sqrt(1:4)
  rcpparma_innerproduct(x)
#> Error in rcpparma_innerproduct(x): could not find function "rcpparma_innerproduct"
  rcpparma_outerproduct(x)
#> Error in rcpparma_outerproduct(x): could not find function "rcpparma_outerproduct"
```
