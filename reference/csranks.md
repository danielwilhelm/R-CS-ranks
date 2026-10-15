# Confidence sets for ranks

Marginal and simultaneous confidence sets for ranks.

## Usage

``` r
csranks(
  x,
  Sigma,
  coverage = 0.95,
  cstype = "two-sided",
  stepdown = TRUE,
  R = 1000,
  simul = TRUE,
  indices = NA,
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

- coverage:

  nominal coverage of the confidence set. Default is 0.95.

- cstype:

  type of confidence set (`two-sided`, `upper`, `lower`). Default is
  `two-sided`.

- stepdown:

  logical; if `TRUE` (default), stepwise procedure is used, otherwise
  single step procedure is used. See Details section for more.

- R:

  number of bootstrap replications. Default is 1000.

- simul:

  logical; if `TRUE` (default), then simultaneous confidence sets are
  computed, which jointly cover all populations indicated by `indices`.
  Otherwise, for each population indicated in `indices` a marginal
  confidence set is computed.

- indices:

  vector of indices of `x` for whose ranks the confidence sets are
  computed. `indices=NA` (default) means computation of confidence sets
  for all populations.

- na.rm:

  logical; if `TRUE`, then `NA`'s are removed from `x` and `Sigma` (if
  any).

- seed:

  seed for bootstrap random variable draws. If set to `NA` (default),
  then seed is not set.

## Value

A `csranks` object, which is a list with three items:

- `L`:

  Lower bounds of the confidence sets for ranks indicated in `indices`

- `rank`:

  Estimated ranks from [`irank`](irank.md) with default parameters

- `U`:

  Upper bounds of the confidence sets.

## Details

Suppose \\j=1,\ldots,p\\ populations (e.g., schools, hospitals,
political parties, countries) are to be ranked according to some measure
\\\theta=(\theta_1,\ldots,\theta_p)\\. We do not observe the true values
\\\theta_1,\ldots,\theta_p\\. Instead, for each population, we have data
from which we have estimated these measures,
\\\hat{\theta}=(\hat{\theta}\_1,\ldots,\hat{\theta}\_p)\\. The values
\\\hat{\theta}\_1,\ldots,\hat{\theta}\_p\\ are estimates of the true
values \\\theta_1,\ldots,\theta_p\\ and thus contain statistical
uncertainty. In consequence, a ranking of the populations by the values
\\\hat{\theta}\_1,\ldots,\hat{\theta}\_p\\ contains statistical
uncertainty and is not necessarily equal to the true ranking of
\\\theta_1,\ldots,\theta_p\\.

The function computes confidence sets for the rank of one, several or
all of the populations (`indices` indicates which of the \\1,\ldots,p\\
populations are of interest). `x` is a vector containing the estimates
\\\hat{\theta}\_1,\ldots,\hat{\theta}\_p\\ and `Sigma` is an estimate of
the covariance matrix of `x`. The method assumes that the estimates are
asymptotically normal and the sample sizes of the datasets are large
enough so that \\\hat{\theta}-\theta\\ is approximately distributed as
\\N(0,\Sigma)\\. The argument `Sigma` should contain an estimate of the
covariance matrix \\\Sigma\\. For instance, if for each population \\j\\
\$\$\sqrt{n_j} (\hat{\theta}\_j-\theta_j) \to_d N(0, \sigma_j^2)\$\$ and
the datasets for each population are drawn independently of each other,
then `Sigma` is a diagonal matrix
\$\$diag(\hat{\sigma}\_1^2/n_1,\ldots,\hat{\sigma}\_p^2/n_p)\$\$
containing estimates of the asymptotic variances divided by the sample
size. More generally, the estimates in `x` may be dependent, but then
`Sigma` must be an estimate of its covariance matrix including
off-diagonal terms.

Marginal confidence sets (`simul=FALSE`) are such that the confidence
set for a population \\j\\ contains the true rank of that population
\\j\\ with probability approximately equal to the nominal coverage
level. Simultaneous confidence sets (`simul=TRUE`) on the other hand are
such that the confidence sets for populations indicated in `indices`
cover the true ranks of all of these populations simultaneously with
probability approximately equal to the nominal coverage level. For
instance, in the PISA example below, a marginal confidence set of a
country \\j\\ covers the true rank of country \\j\\ with probability
approximately equal to 0.95. A simultaneous confidence set for all
countries covers the true ranks of all countries simultaneously with
probability approximately equal to 0.95.

The function implements the procedures developed and described in more
detail in Mogstad, Romano, Shaikh, and Wilhelm (2023). The procedure is
based on on testing a large family of hypotheses for pairwise
comparisons. Stepwise methods can be used to improve the power of the
procedure by, potentially, rejecting more hypotheses without violating
the desired coverage property of the resulting confidence set. These are
employed when `stepdown=TRUE`. From a practical point of view,
`stepdown=TRUE` is computationally more demanding, but often results in
tighter confidence sets.

The procedure uses a parametric bootstrap procedure based on the above
approximate multivariate normal distribution.

## References

Mogstad, Romano, Shaikh, and Wilhelm (2023), "Inference for Ranks with
Applications to Mobility across Neighborhoods and Academic Achievements
across Countries", forthcoming at Review of Economic Studies [cemmap
working paper](http://dwilhelm.userweb.mwn.de/papers/cwp0323.pdf)
[doi:10.1093/restud/rdad006](https://doi.org/10.1093/restud/rdad006)

## Examples

``` r
# simple simulated example:
n <- 100
p <- 10
X <- matrix(rep(1:p, n) / p, ncol = p, byrow = TRUE) + matrix(rnorm(n * p), 100, 10)
thetahat <- colMeans(X)
Sigmahat <- cov(X) / n
csranks(thetahat, Sigmahat)
#> $L
#>  [1] 4 4 4 2 4 2 2 1 1 1
#> 
#> $rank
#>  [1] 10  8  9  4  7  5  6  2  3  1
#> 
#> $U
#>  [1] 10 10 10 10 10 10 10  6  6  3
#> 
#> attr(,"class")
#> [1] "csranks"

# PISA example:
data(pisa2018)
math_score <- pisa2018$math_score
math_se <- pisa2018$math_se
math_cov_mat <- diag(math_se^2)

# marginal confidence set for each country:
csranks(math_score, math_cov_mat, simul = FALSE)
#> $L
#>  [1] 15 10  5  4 35 37 36 10  5  1  5 12 10 32 25 12 10 32 19  1  1 12 25 25 35
#> [26]  1 12 10  2 12 20  5 25  7  2 32  8 25
#> 
#> $rank
#>  [1] 24 18 10  7 35 38 37 17  8  3 11 20 15 34 30 21 16 32 25  1  2 19 29 27 36
#> [26]  4 22 14  5 23 26  9 28 12  6 33 13 31
#> 
#> $U
#>  [1] 26 24 18 13 36 38 38 24 13  6 18 26 24 34 31 26 24 34 31  4  6 25 31 31 37
#> [26]  7 26 23 11 26 31 14 31 23 11 34 23 31
#> 
#> attr(,"class")
#> [1] "csranks"

# simultaneous confidence set for all countries:
csranks(math_score, math_cov_mat, simul = TRUE)
#> $L
#>  [1] 12  7  4  3 35 37 35  7  4  1  4 12  7 32 23 12  8 31 15  1  1 12 23 23 35
#> [26]  1 12  7  1 12 18  4 24  5  1 32  6 23
#> 
#> $rank
#>  [1] 24 18 10  7 35 38 37 17  8  3 11 20 15 34 30 21 16 32 25  1  2 19 29 27 36
#> [26]  4 22 14  5 23 26  9 28 12  6 33 13 31
#> 
#> $U
#>  [1] 30 26 18 17 37 38 38 25 18  6 18 26 25 34 31 26 25 34 31  6  7 26 31 31 37
#> [26] 11 26 24 12 31 31 18 31 24 13 34 24 32
#> 
#> attr(,"class")
#> [1] "csranks"
```
