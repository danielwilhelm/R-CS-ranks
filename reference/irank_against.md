# Compute integer ranks in another reference vector

The method [`irank`](irank.md) compares ranks using the same vector as
reference. `irank_against` returns integer ranks, that values from `x`
would assume if (individually) inserted into `v`. `frank_against` acts
analogously, returning fractional ranks.

## Usage

``` r
irank_against(
  x,
  v,
  omega = 0,
  increasing = FALSE,
  na.rm = FALSE,
  weights = NULL
)

frank_against(
  x,
  v,
  omega = 0,
  increasing = FALSE,
  na.rm = FALSE,
  weights = NULL
)
```

## Arguments

- x:

  numeric query vector.

- v:

  numeric reference vector.

- omega:

  numeric value in \[0,1\], defining how ties in `x` (if any) are
  handled; default is `0`. See Details.

- increasing:

  logical; if `FALSE` (default), then large elements in `x` receive a
  small rank. Otherwise, large elements in `x` receive a large rank.

- na.rm:

  logical; if `TRUE`, then `NA`'s are removed from `x`. Default:
  `FALSE`.

- weights:

  an optional vector of weights to be used in the ranking process.
  Should be `NULL` or a numeric vector.

## Value

Numeric vector of the same length as `x` containing the integer (for
`irank_against`) or fractional (for `frank_against`) ranks.

## Details

It's useful to think about `frank_against(x,v)` as a generalization of
Empirical Cumulative Distribution Function, created for `v` and
evaluated for points in `x`. `frank_agaist(x,v,increasing=TRUE,omega=1)`
is identical to `ecdf(v)(x)`.

`increasing` switches the inequality sign in ECDF definition from
\\F_V(t) = \hat P(V \<= t)\\ to \\\hat P(V \>= t)\\.

`omega=0` introduces the strict inequality (\\\hat P(V \< t)\\ instead
of \\\hat P(V \<= t)\\). Any `omega` in between is a weighted average of
the cases `omega=1` and `omega=0`.

Finally, `irank_against` is equal to `frank_against` multiplied by the
`length(v)`.

This particular choice of default parameters was made for compatibility
with default parameters of `irank` and `frank`. `irank(x)` is always
equal to `irank_against(x,x)` and `frank(x)` is always equal to
`frank_against(x,x)`.

## See also

[`irank()`](irank.md), [`ecdf()`](https://rdrr.io/r/stats/ecdf.html)

## Examples

``` r
irank_against(1:10, c(4, 4, 4, 3, 1, 10, 7, 7))
#>  [1] 8 8 7 4 4 4 2 2 2 1
```
