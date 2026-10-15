# Confidence sets for the tau-best

Computation of confidence sets for the identities of populations among
the tau best.

## Usage

``` r
cstaubest(
  x,
  Sigma,
  tau = 2,
  coverage = 0.95,
  stepdown = TRUE,
  R = 1000,
  na.rm = FALSE,
  seed = NA
)

cstauworst(
  x,
  Sigma,
  tau = 2,
  coverage = 0.95,
  stepdown = TRUE,
  R = 1000,
  na.rm = FALSE,
  seed = NA
)
```

## Arguments

- x:

  vector of estimates containing estimated features by which the
  populations are to be ranked.

- Sigma:

  estimated covariance matrix of `x`.

- tau:

  the confidence set contains indicators for the elements in `x` whose
  rank is less than or equal to `tau`.

- coverage:

  nominal coverage of the confidence set. Default is 0.95.

- stepdown:

  logical; if `TRUE` (default), stepwise procedure is used, otherwise
  single step procedure is used. See Details section for more.

- R:

  number of bootstrap replications. Default is 1000.

- na.rm:

  logical; if `TRUE`, then `NA`'s are removed from `x` and `Sigma` (if
  any).

- seed:

  seed for bootstrap random variable draws. If set to `NA` (default),
  then seed is not set.

## Value

logical vector indicating which of the elements of `x` are in the
confidence set for the tau-best.

## Functions

- `cstauworst()`: Confidence sets for the tau-worst

  Equivalent to calling `cstaubest` with `-x`.

## Details

The function computes a confidence set containing indicators for the
elements in `x` whose rank is less than or equal to `tau` with
probability approximately equal to the nominal coverage (`coverage`).

The function implements the projection confidence set for the tau-best
developed and described in more detail in Mogstad, Romano, Shaikh, and
Wilhelm (2023).

## References

Mogstad, Romano, Shaikh, and Wilhelm (2023), "Inference for Ranks with
Applications to Mobility across Neighborhoods and Academic Achievements
across Countries", forthcoming at Review of Economic Studies [cemmap
working paper](http://dwilhelm.userweb.mwn.de/papers/cwp0323.pdf),
[doi:10.1093/restud/rdad006](https://doi.org/10.1093/restud/rdad006)

## Examples

``` r
# simple simulated example:
n <- 100
p <- 10
X <- matrix(rep(1:p, n) / p, ncol = p, byrow = TRUE) + matrix(rnorm(n * p), 100, 10)
thetahat <- colMeans(X)
Sigmahat <- cov(X) / n

# confidence set for the populations that may be among the top-3
# (with probability approximately 0.95):
cstaubest(thetahat, Sigmahat, tau = 3)
#>  [1] FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE

# confidence set for the populations that may be among the bottom-3
# (with probability approximately 0.95):
cstauworst(thetahat, Sigmahat, tau = 3)
#>  [1]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE
```
