
<!-- README.md is generated from README.Rmd. Please edit that file -->

# csranks

<!-- badges: start -->

[![R-CMD-check](https://github.com/danielwilhelm/R-CS-ranks/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/danielwilhelm/R-CS-ranks/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The `R` package `csranks` provides statistical tools for estimation and
inference involving ranks (the position in a ranking). Two central
functions are `csranks` for confidence sets for ranks and `lmranks` for
regressions involving ranks, e.g. rank-rank regressions that are popular
in applied work in economics.

The functions are based on recent work developing these procedures and
their theoretical properties. The confidence sets for ranks are based on
[Mogstad, Romano, Shaikh, and Wilhelm
(2023)](https://doi.org/10.1093/restud/rdad006) and [Bazylik, Mogstad,
Romano, Shaikh, and Wilhelm
(2022)](https://dwilhelm.userweb.mwn.de/papers/cwp4021.pdf). The
inference methods for regressions involving ranks are developed in
[Chetverikov and Wilhelm (2023)](https://arxiv.org/pdf/2310.15512).

## Installation

You can install the released version of `csranks` from CRAN with:

``` r
install.packages("csranks")
```

You can install the development version of `csranks` from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("danielwilhelm/R-CS-ranks")
```
