# Compute ranks

Compute integer of fractional ranks with flexible handling of ties.

## Usage

``` r
irank(x, omega = 0, increasing = FALSE, na.rm = FALSE, weights = NULL)

frank(x, omega = 0, increasing = FALSE, na.rm = FALSE, weights = NULL)
```

## Arguments

- x:

  vector of values to be ranked

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
`irank`) or fractional (for `frank`) ranks.

## Details

`irank` implements all possible definitions of ranks of the values in
`x`. Different definitions of the ranks are chosen through combinations
of the two arguments `omega` and `increasing`. Suppose `x` is of length
\\p\\. If `increasing=TRUE`, then the largest value in `x` receives the
rank \\p\\ and the smallest the rank \\1\\. If `increasing=FALSE`, then
the largest value in `x` receives the rank \\1\\ and the smallest the
rank \\p\\.

The value of `omega` indicates how ties are handled. If there are no
ties in `x`, then the value of `omega` does not affect the ranks and the
only choice to be made is whether the ranks should be increasing or
decreasing with the values in `x`. When there are ties in `x`, however,
then there are infinitely many possible ranks that can be assigned to a
tied value.

When `increasing=TRUE`, then `omega=0` leads to the smallest possible
and `omega=1` to the largest possible rank of a tied value. Values of
`omega` between 0 and 1 lead to values of the rank between the largest
and smallest.

Gotcha: the definition of weighted ranks drops a constant (1-`omega`).
This means that `irank` (and other functions from the family) called
with equal weights will have different outcome from `irank` called with
`weights=NULL`. The values will be shifted by a constant 1-`omega`.

`frank` takes the ranking returned by `irank` and divides the result by
`length(x)`. The result is a ranking with ranks in the interval \[0,1\].
An important special case occurs for `increasing=TRUE` and `omega=1`: in
this case, the rank of the value `x[j]` is equal to the empirical cdf of
`x` evaluated at `x[j]`.

## Examples

``` r
# simple example without ties:
x <- c(3, 8, -4, 10, 2)
irank(x, increasing = TRUE)
#> [1] 3 4 1 5 2
irank(x, increasing = FALSE)
#> [1] 3 2 5 1 4

# since there are no ties, the value of omega has no impact:
irank(x, increasing = TRUE, omega = 0)
#> [1] 3 4 1 5 2
irank(x, increasing = TRUE, omega = 0.5)
#> [1] 3 4 1 5 2
irank(x, increasing = TRUE, omega = 1)
#> [1] 3 4 1 5 2

# simple example with ties:
x <- c(3, 4, 7, 7, 10, 11, 15, 15, 15, 15)
irank(x, increasing = TRUE, omega = 0) # smallest possible ranks
#>  [1] 1 2 3 3 5 6 7 7 7 7
irank(x, increasing = TRUE, omega = 0.5) # mid-ranks
#>  [1] 1.0 2.0 3.5 3.5 5.0 6.0 8.5 8.5 8.5 8.5
irank(x, increasing = TRUE, omega = 1) # largest possible ranks
#>  [1]  1  2  4  4  5  6 10 10 10 10

# simple example of fractional ranks without ties:
x <- c(3, 8, -4, 10, 2)
frank(x, increasing = TRUE)
#> [1] 0.6 0.8 0.2 1.0 0.4
frank(x, increasing = FALSE)
#> [1] 0.6 0.4 1.0 0.2 0.8
```
