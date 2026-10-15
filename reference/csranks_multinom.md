# Confidence sets for ranks based on multinomial data

Marginal and simultaneous confidence sets for ranks of categories, where
categories are ranked by the probabilities of being chosen.

## Usage

``` r
csranks_multinom(
  x,
  coverage = 0.95,
  cstype = "two-sided",
  simul = TRUE,
  multcorr = "Holm",
  indices = NA,
  na.rm = FALSE
)
```

## Arguments

- x:

  vector of counts indicating how often each category was chosen.

- coverage:

  nominal coverage of the confidence set. Default is 0.95.

- cstype:

  type of confidence set (`two-sided`, `upper`, `lower`). Default is
  `two-sided`.

- simul:

  logical; if `TRUE` (default), then simultaneous confidence sets are
  computed, which jointly cover all populations indicated by `indices`.
  Otherwise, for each population indicated in `indices` a marginal
  confidence set is computed.

- multcorr:

  multiplicity correction to be used: `Holm` (default) or `Bonferroni`.
  See Details section for more.

- indices:

  vector of indices of `x` for whose ranks the confidence sets are
  computed. `indices=NA` (default) means computation of confidence sets
  for all populations.

- na.rm:

  logical; if `TRUE`, then `NA`'s are removed from `x` and `Sigma` (if
  any).

## Value

A `csranks` object, which is a list with three items:

- `L`:

  Lower bounds of the confidence sets for ranks indicated in `indices`

- `rank`:

  Estimated ranks from [`irank`](irank.md) with default parameters

- `U`:

  Upper bounds of the confidence sets.

## Details

This function computes confidence sets for ranks similarly as
[`csranks`](csranks.md), but it is tailored to the special case of
multinomial data. Suppose there are \\p\\ populations (for the case of
multinomial data, we will refer to them as "categories") such as
political parties, for example, that one wants to rank by the
probabilities of them being chosen. For political parties, this would
correspond to the share of votes each party obtains. Here, the
underlying data are multinomial: each observation corresponds to a
choice among the \\p\\ categories. The vector `x` contains the counts of
how often each category was chosen in the data.

In this setting, `link{csranks}` could be applied to compute confidence
sets for the ranks of each category, but instead this function
implements a different method proposed by Bazylik, Mogstad, Romano,
Shaikh, and Wilhelm (2023), which exploits the multinomial structure of
the problem and yields confidence sets for the ranks that are valid in
finite samples (whereas [`csranks`](csranks.md) produces confidence sets
that are valid only asymptotically).

The procedure involves testing multiple hypotheses. The
`\code{multcorr}` indicates a method for multiplicity correction. See
the paper for details.

## References

Bazylik, Mogstad, Romano, Shaikh, and Wilhelm. "Finite-and large-sample
inference for ranks using multinomial data with an application to
ranking political parties". [cemmap working
paper](http://dwilhelm.userweb.mwn.de/papers/cwp4021.pdf)

## Examples

``` r
x <- c(rmultinom(1, 1000, 1:10))
csranks_multinom(x)
#> $L
#>  [1] 9 8 6 5 5 2 2 2 1 1
#> 
#> $rank
#>  [1] 10  9  8  7  6  5  4  3  2  1
#> 
#> $U
#>  [1] 10 10  9  8  8  7  5  5  5  2
#> 
#> attr(,"class")
#> [1] "csranks"
```
