# csranks (development version)
* Fixed S3 inheritance for the `proj` method of `lmranks`. Added it to `lmranks` documentation page.
* Changed column names of `coefficients` in `summary.lmranks` from `t value` to `z value` and from `Pr(<|t|)` to `Pr(<|z|)`.
* Removed line about residual degrees of freedom from `print.summary.lmranks`.
* Example in `csranks` now calls `data(pisa2018)`.
* Clarified documentation of `csranks` function family.
* Improved documentation of `irank`, `frank`, etc
* Fixed typo in an error message
* Updated tests

# csranks 1.2.2
* Changed license to GPLv3

# csranks 1.2.1

* Added new datasets, `pisa2018` and `pisa2022`
* Updated CITATION and references
* Bugfix of `plot.csranks` with non default `indices`

# csranks 1.2.0
* Public release with contents of 1.1 (Rank-rank linear regression)
* Added a new mode for regression with clustered (grouped) data
* Bugfixes, doc, new vignette

# csranks 1.1.1
* Optimized `vcov.lmranks` method

# csranks 1.1.0

* Added `lmranks` function for linear modelling of ranks using single rank
covariate and possibly other, "usual" covariates
* Implemented methods `print`, `summary`, `vcov`, `confint`, `predict` for `lmranks` output
* Disabled a number of methods defined for `ļm` (like `sigma`, `AIC` or `influence`)

# csranks 1.0.0

* Release!
* `irank` now raises error if NAs present and `na.rm`=FALSE

# csranks 0.5.0

* `csranks` now produces `csranks` object. List as before, but with new `rank` element
* Added `plot.ranking` function

# csranks 0.4.1

* `V` argument (`sd` before 0.4) renamed to `Sigma`
* Adding `irank` and `frank` functions

# csranks 0.4.0

* `sd` argument renamed to `V`; now accepts ONLY covariance matrix.
* `*_marg` and `*_simul` methods removed; use `*` with `simul` option.

# csranks 0.3.0

* Added possibility of features to be correlated across populations. 
`sd` argument now can accept covariance matrix.

# csranks 0.2.0

* Added implementation of confidence sets for ranks based on multinomial data.
* Added implementation of confidence sets for tau-best and tau-worst.

# csranks 0.1.0

* This is the first release of csranks.
