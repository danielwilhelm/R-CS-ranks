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
