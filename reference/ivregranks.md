# Instrumental-Variable Regression by 2SLS Estimation Involving Ranks

Fit instrumental-variable regression involving ranks by two-stage least
squares (2SLS). This is equivalent to direct instrumental-variables
estimation when the number of instruments is equal to the number of
regressors.

## Usage

``` r
ivregranks(
  formula,
  data,
  subset,
  na.action,
  weights,
  offset,
  contrasts = NULL,
  model = TRUE,
  y = TRUE,
  x = FALSE,
  method = "OLS",
  omega = 1,
  ...
)

# S3 method for class 'ivregranks'
plot(x, which = 1, ...)

# S3 method for class 'ivregranks'
predict(object, newdata, ...)

# S3 method for class 'ivregranks'
summary(object, vcov. = NULL, df = NULL, diagnostics = NULL, ...)

# S3 method for class 'summary.ivregranks'
print(x, ...)

# S3 method for class 'ivregranks'
confint(
  object,
  parm,
  level = 0.95,
  component = c("stage2", "stage1"),
  complete = TRUE,
  vcov. = NULL,
  df = NULL,
  ...
)

# S3 method for class 'ivregranks'
vcov(object, component = c("stage2", "stage1"), complete = TRUE, ...)
```

## Arguments

- formula:

  formula specification(s) of the regression relationship and the
  instruments. `formula` has three parts as in
  `r(y) ~ x1 + r(x2) | r(z1) + z2 + z3`.

- data:

  an optional data frame containing the variables in the model. By
  default the variables are taken from the environment of the `formula`.

- subset:

  currently not supported.

- na.action:

  currently not supported.

- weights:

  currently not supported.

- offset:

  an optional offset that can be used to specify an a priori known
  component to be included during fitting.

- contrasts:

  an optional list. See the `contrasts.arg` of
  [`model.matrix.default`](https://rdrr.io/r/stats/model.matrix.html).

- model, x, y:

  logicals. If `TRUE` the corresponding components of the fit (the model
  frame, the model matrices, the response) are returned. These
  components are necessary for computing regression diagnostics.

- method:

  the method used to fit the stage 1 and 2 regression: `"OLS"` for
  traditional 2SLS regression (the default and only option).

- omega:

  real number in the interval \[0,1\] defining how ties are handled (if
  there are any).

- ...:

  further arguments passed to
  [`ivreg.fit`](https://zeileis.github.io/ivreg/reference/ivreg.fit.html).

- which:

  As in
  [`plot.ivreg`](https://zeileis.github.io/ivreg/reference/ivregDiagnostics.html).
  Currently only no. 1 is available.

- object:

  An object of class `"ivregranks"`.

- newdata:

  An optional data frame in which to look for variables with which to
  predict. If omitted, the fitted values are used.

- vcov.:

  Optionally either a coefficient covariance matrix or a function to
  compute such a covariance matrix from fitted `ivreg` model objects. If
  `NULL` (the default) the standard covariance matrix (based on the
  information matrix) is used. Alternatively, covariance matrices (e.g.,
  clustered and/or heteroscedasticity-consistent) can be plugged in to
  adjust Wald tests or confidence intervals etc. In `summary`, if
  `diagnostics = TRUE`, `vcov.` must be a function (not a matrix)
  because the alternative covariances are also needed for certain
  auxiliary models in the diagnostic tests. If `vcov.` is a function,
  the `...` argument can be used to pass on further arguments to this
  function.

- df:

  For `summary`, optional residual degrees of freedom to use in
  computing model summary.

- diagnostics:

  currently not supported.

- parm:

  parameters for which confidence intervals are to be computed; a vector
  or numbers or names; the default is all parameters.

- level:

  confidence level; the default is `0.95`.

- component:

  Character indicating `"stage2"` or `"stage1"`.

- complete:

  logical indicating if the full variance-covariance matrix should be
  returned also in case of an over-determined system where some
  coefficients are undefined and `coef(.)` contains NAs correspondingly.
  When `complete = TRUE`, [`vcov()`](https://rdrr.io/r/stats/vcov.html)
  is compatible with [`coef()`](https://rdrr.io/r/stats/coef.html) also
  in this singular case.

## Value

`ivregranks` returns an object of class `"ivregranks"` that inherits as
much as possible from class
[`ivreg`](https://zeileis.github.io/ivreg/reference/ivreg.html), with
the following additional components:

- rank_terms_indices:

  an integer vector with indices of entries of `terms.labels` attribute
  of `terms(formula)` for the outcome equation which correspond to
  ranked regressors.

- rank_instruments_indices:

  an integer vector with indices of entries of the ranked instrumental
  variables.

- ranked_response:

  a logical entry.

- omega:

  an entry corresponding to the `omega` argument.

## Details

Regressors and instruments for `ivregranks` are most easily specified in
a formula with two parts on the right-hand side, e.g.,
`r(y) ~ x1 + r(x2) | r(z1) + z2 + z3`, where `x1` and `r(x2)` are the
explanatory variables and `r(z1)`, `z2`, and `z3` are the instrumental
variables. Note that exogenous regressors have to be included as
instruments for themselves.

For example, if there is one exogenous regressor `ex` and one endogenous
regressor `r(en)` with instrument `r(in)`, the appropriate formula would
be `r(y) ~ r(en) + ex | r(in) + ex`. Alternatively, a formula with three
parts on the right-hand side can also be used:
`r(y) ~ ex | r(en) | r(in)`. The latter is typically more convenient, if
there is a large number of exogenous regressors.

Moreover, two further equivalent specification strategies are possible
that are typically less convenient compared to the strategies above. One
option is to use an update formula with a `.` in the second part of the
formula is used: `r(y) ~ r(en) + ex | . - r(en) + r(in)`. Another option
is to use a separate formula for the instruments (only for backward
compatibility with earlier versions):
`formula = r(y) ~ r(en) + ex, instruments = ~ r(in) + ex`.

Internally, all specifications are converted to the version with two
parts on the right-hand side.

## Methods (by generic)

- `plot(ivregranks)`: Plot diagnostics for an `ivregranks` object

  Displays plots useful for assessing quality of model fit. Currently,
  only one plot is available, which plots fitted values against
  residuals (for homoscedacity check).

- `predict(ivregranks)`: Predict method for IV Model for Ranks Fits

- `summary(ivregranks)`: Summary and Inference Methods for
  `"ivregranks"` Objects

- `vcov(ivregranks)`: Calculate Variance-Covariance Matrix for a Fitted
  `ivregranks` object

  Returns the variance-covariance matrix of the regression coefficients
  (main parameters) of a fitted `ivregranks` object. Its result is
  theoretically valid and asymptotically consistent, in contrast to
  naively running `vcov(ivreg(...))`.

## References

Chetverikov and Wilhelm (2023), "Inference for Rank-Rank Regressions".
[arXiv preprint arXiv:2310.15512](https://arxiv.org/pdf/2310.15512)

## See also

[`ivreg.fit`](https://zeileis.github.io/ivreg/reference/ivreg.fit.html),
[`lmranks`](lmranks.md)

Generic functions [`coef`](https://rdrr.io/r/stats/coef.html),
[`residuals`](https://rdrr.io/r/stats/residuals.html),
[`fitted`](https://rdrr.io/r/stats/fitted.values.html),
[`model.frame`](https://rdrr.io/r/stats/model.frame.html),
[`model.matrix`](https://rdrr.io/r/stats/model.matrix.html),
[`update`](https://rdrr.io/r/stats/update.html) .

## Examples

``` r
# rank-rank regression:
Z <- rnorm(500)
X <- Z + rnorm(500)
Y <- X + rnorm(500)
rrfit <- ivregranks(r(Y) ~ r(X) | r(Z))
summary(rrfit)
#> 
#> Call:
#> ivregranks(formula = r(Y) ~ r(X) | r(Z))
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -0.535939 -0.124013  0.002403  0.114031  0.556328 
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.10579    0.01658   6.379 1.78e-10 ***
#> r(X)         0.78884    0.03310  23.831  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 

# naive version of the rank-rank regression:
RZ <- frank(Z, increasing = TRUE, omega = 1)
RX <- frank(X, increasing = TRUE, omega = 1)
RY <- frank(Y, increasing = TRUE, omega = 1)
fit <- ivreg::ivreg(RY ~ RX | RZ)
summary(fit)
#> 
#> Call:
#> ivreg::ivreg(formula = RY ~ RX | RZ)
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -0.535939 -0.124013  0.002403  0.114031  0.556328 
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  0.10579    0.02025   5.225 2.57e-07 ***
#> RX           0.78884    0.03723  21.190  < 2e-16 ***
#> 
#> Diagnostic tests:
#>                  df1 df2 statistic p-value    
#> Weak instruments   1 498   580.619  <2e-16 ***
#> Wu-Hausman         1 497     0.024   0.876    
#> Sargan             0  NA        NA      NA    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual standard error: 0.1763 on 498 degrees of freedom
#> Multiple R-Squared: 0.6285,  Adjusted R-squared: 0.6278 
#> Wald test:   449 on 1 and 498 DF,  p-value: < 2.2e-16 
#> 
# the coefficient estimates are the same as in the ivregranks function, but
# the standard errors, t-values, p-values are incorrect.

# support for `data` argument:
ivr <- ivregranks(r(mpg) ~ r(hp) + cyl | r(disp) + cyl, data = mtcars)
summary(ivr)
#> 
#> Call:
#> ivregranks(formula = r(mpg) ~ r(hp) + cyl | r(disp) + cyl, data = mtcars)
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -1.49398 -0.41639 -0.06979  0.44197  2.02883 
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)
#> (Intercept)  -1.4537     9.5927  -0.152    0.880
#> r(hp)        -7.1939    23.4357  -0.307    0.759
#> cyl           0.9336     3.5509   0.263    0.793
#> 
```
