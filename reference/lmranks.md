# Regressions Involving Ranks

Estimation and inference for regressions involving ranks, i.e.
regressions in which the dependent and/or the independent variable has
been transformed into ranks before running the regression.

## Usage

``` r
lmranks(
  formula,
  data,
  subset,
  weights,
  na.action = stats::na.fail,
  method = "qr",
  model = TRUE,
  x = FALSE,
  qr = TRUE,
  y = FALSE,
  singular.ok = TRUE,
  contrasts = NULL,
  offset = offset,
  omega = 1,
  ...
)

# S3 method for class 'lmranks'
plot(x, which = 1, ...)

# S3 method for class 'lmranks'
proj(object, onedf = FALSE, ...)

# S3 method for class 'lmranks'
predict(object, newdata, ...)

# S3 method for class 'lmranks'
summary(object, correlation = FALSE, symbolic.cor = FALSE, ...)

# S3 method for class 'lmranks'
vcov(object, complete = TRUE, ...)
```

## Arguments

- formula:

  An object of class
  "[`formula`](https://rdrr.io/r/stats/formula.html)": a symbolic
  description of the model to be fitted. Exactly like the formula for
  linear model except that variables to be ranked can be indicated by
  `r()`. See Details and Examples below.

- data:

  an optional data frame, list or environment (or object coercible by
  [`as.data.frame`](https://rdrr.io/r/base/as.data.frame.html) to a data
  frame) containing the variables in the model. If not found in `data`,
  the variables are taken from `environment(formula)`, typically the
  environment from which `lm` is called.

- subset:

  currently not supported.

- weights:

  an optional vector of weights to be used in the fitting process.
  Should be `NULL` or a numeric vector. If non-NULL, weighted least
  squares is used with weights `weights` (that is, minimizing
  `sum(w*e^2)`); otherwise ordinary least squares is used.

- na.action:

  currently not supported. User is expected to handle NA values prior to
  the use of this function.

- method:

  the method to be used; for fitting, currently only `method = "qr"` is
  supported; `method = "model.frame"` returns the model frame (the same
  as with `model = TRUE`, see below).

- model, y, qr:

  logicals. If TRUE the corresponding components of the fit (the model
  frame, the response, the QR decomposition) are returned.

- x:

  - For `lmranks`: Logical. Should model matrix be returned?

  - For `plot` method: An `lmranks` object.

- singular.ok:

  logical. If `FALSE` (the default in S but not in R) a singular fit is
  an error.

- contrasts:

  an optional list. See the `contrasts.arg` of
  [`model.matrix.default`](https://rdrr.io/r/stats/model.matrix.html).

- offset:

  this can be used to specify an *a priori* known component to be
  included in the linear predictor during fitting. This should be `NULL`
  or a numeric vector or matrix of extents matching those of the
  response. One or more [`offset`](https://rdrr.io/r/stats/offset.html)
  terms can be included in the formula instead or as well, and if more
  than one are specified their sum is used. See
  [`model.offset`](https://rdrr.io/r/stats/model.extract.html).

- omega:

  real number in the interval \[0,1\] defining how ties are handled (if
  there are any). The value of `omega` is passed to [`frank`](irank.md)
  for computation of ranks. The default is 1 so that the rank of a
  realized value is defined as the empirical cdf evaluated at that
  realized value. See Details below.

- ...:

  For [`lm()`](https://rdrr.io/r/stats/lm.html): additional arguments to
  be passed to the low level regression fitting functions (see below).

- which:

  As in [`plot.lm`](https://rdrr.io/r/stats/plot.lm.html). Currently
  only no.1 is available.

- object:

  A `lmranks` object.

- onedf:

  A logical flag. If `TRUE`, a projection is returned for all the
  columns of the model matrix. If `FALSE`, the single-column projections
  are collapsed by terms of the model (as represented in the analysis of
  variance table).

- newdata:

  An optional data frame in which to look for variables with which to
  predict. If omitted, the fitted values are used.

- correlation:

  logical; if `TRUE`, the correlation matrix of the estimated parameters
  is returned and printed.

- symbolic.cor:

  logical. If `TRUE`, print the correlations in a symbolic form (see
  [`symnum`](https://rdrr.io/r/stats/symnum.html)) rather than as
  numbers.

- complete:

  logical indicating if the full variance-covariance matrix should be
  returned also in case of an over-determined system where some
  coefficients are undefined and `coef(.)` contains NAs correspondingly.
  When `complete = TRUE`, [`vcov()`](https://rdrr.io/r/stats/vcov.html)
  is compatible with [`coef()`](https://rdrr.io/r/stats/coef.html) also
  in this singular case.

## Value

An object of class `lmranks`, inheriting (as much as possible) from
class `lm`.

Additionally, it has an `omega` entry, corresponding to the `omega`
argument, a `ranked_response` logical entry, and a
`rank_terms_indices` - an integer vector with indices of entries of
`terms.labels` attribute of `terms(formula)`, which correspond to ranked
regressors.

## Details

This function performs estimation and inference for regressions
involving ranks. Suppose there is a dependent variable \\Y_i\\ and
independent variables \\X_i\\ and \\W_i\\, where \\X_i\\ is a scalar and
\\W_i\\ a vector (possibly including a constant). Instead of running a
linear regression of \\Y_i\\ on \\X_i\\ and \\W_i\\, we want to first
transform \\Y_i\\ and/or \\X_i\\ into ranks. Denote by \\R_i^Y\\ the
rank of \\Y_i\\ and \\R_i^X\\ the rank of \\X_i\\. Then, a **rank-rank
regression**, \$\$R_i^Y = \rho R_i^X + W_i'\beta + \varepsilon_i,\$\$ is
run using the formula `r(Y)~r(X)+W`. Similarly, a **regression of the
raw dependent variable on the ranked regressor**, \$\$Y_i = \rho R_i^X +
W_i'\beta + \varepsilon_i,\$\$ can be implemented by the formula
`Y~r(X)+W`, and a **regression of the ranked dependent variable on the
raw regressors**, \$\$R^Y_i = W_i'\beta + \varepsilon_i,\$\$ can be
implemented by the formula `r(Y)~W`.

The function works, in many ways, just like `lm` for linear regressions.
Apart from some smaller details, there are two important differences:
first, in `lmranks`, the mark `r()` can be used in formulas to indicate
variables to be ranked before running the regression and, second,
subsequent use of `summary` produces a summary table with the correct
standard errors, t-values and p-values (while those of the `lm` are not
correct for regressions involving ranks). See Chetverikov and Wilhelm
(2023) for more details.

Many other aspects of the function are similar to `lm`. For instance,
`.` in a formula means 'all columns not otherwise in the formula' just
as in `lm`. An intercept is included by default. In a model specified as
`r(Y)~r(X)+.`, both `r(X)` and `X` will be included in the model - as it
would have been in `lm` and, say,
[`log()`](https://rdrr.io/r/base/Log.html) instead of `r()`. One can
exclude `X` with a `-`, i.e. `r(Y)~r(X)+.-X`. See
[`formula`](https://rdrr.io/r/stats/formula.html) for more about model
specification.

The `r()` is a private alias for [`frank`](irank.md). The `increasing`
argument, provided at individual regressor level, specifies whether the
ranks should increase or decrease as regressor values increase. The
`omega` argument of [`frank`](irank.md), provided at `lmranks` function
level, specifies how ties in variables are to be handled and can be
supplied as argument in `lmranks`. For more details, see
[`frank`](irank.md). By default `increasing` is set to `TRUE` and
`omega` is set equal to `1`, which means `r()` computes ranks by
transforming a variable through its empirical cdf.

Many functions defined for `lm` also work correctly with `lmranks`.
These include [`coef`](https://rdrr.io/r/stats/coef.html),
[`model.frame`](https://rdrr.io/r/stats/model.frame.html),
[`model.matrix`](https://rdrr.io/r/stats/model.matrix.html),
[`resid`](https://rdrr.io/r/stats/residuals.html),
[`update`](https://rdrr.io/r/stats/update.html) and others. On the other
hand, some would return incorrect results if they treated `lmranks`
output in the same way as `lm`'s. The central contribution of this
package are `vcov`, `summary` and `confint` implementations using the
correct asymptotic theory for regressions involving ranks.

See the [`lm`](https://rdrr.io/r/stats/lm.html) documentation for more.

## Methods (by generic)

- `plot(lmranks)`: Plot diagnostics for an `lmranks` object

  Displays plots useful for assessing quality of model fit. Currently,
  only one plot is available, which plots fitted values against
  residuals (for homoscedacity check).

- `proj(lmranks)`: Projections of the data onto terms of rank-rank
  regression model

- `predict(lmranks)`: Predict method for Linear Model for Ranks Fits

- `summary(lmranks)`: Summarizing fits of rank-rank regressions

- `vcov(lmranks)`: Calculate Variance-Covariance Matrix for a Fitted
  `lmranks` object

  Returns the variance-covariance matrix of the regression coefficients
  (main parameters) of a fitted `lmranks` object. Its result is
  theoretically valid and asymptotically consistent, in contrast to
  naively running `vcov(lm(...))`.

## Rank-rank regressions with clusters

Sometimes, the data is divided into clusters (groups) and one is
interested in running rank-rank regressions separately within each
cluster, where the ranks are not computed within each cluster, but using
all observations pooled across all clusters. Specifically, let
\\G_i=1,\ldots,n_G\\ denote a variable that indicates the cluster to
which the i-th observation belongs. Then, the regression model of
interest is \$\$R_i^Y = \sum\_{g=1}^{n_G} 1\\G_i=g\\(\rho_g R_i^X +
W_i'\beta_g) + \varepsilon_i,\$\$ where \\\rho_g\\ and \\\beta_g\\ are
now cluster-specific coefficients, but the ranks \\R_i^Y\\ and \\R_i^X\\
are computed as ranks among all observations \\Y_i\\ and \\X_i\\,
respectively. That means the rank of an observation is not computed
among the other observations in the same cluster, but rather among all
available observations across all clusters.

This type of regression is implemented in the `lmranks` function using
interaction notation: `r(Y)~(r(X)+W):G`. Here, the variable G **must**
be a [`factor`](https://rdrr.io/r/base/factor.html).

Since the theory for clustered regression mixing grouped and ungrouped
(in)dependent variables is not yet developed, such a model will raise an
error. Also, by default the function includes a cluster-specific
intercept, i.e. `r(Y)~(r(X)+W):G` is internally interpreted as
`r(Y)~(r(X)+W):G+G-1`.

[`contrasts`](https://rdrr.io/r/stats/contrasts.html) of `G` must be of
`contr.treatment` kind, which is the default.

## Warning

As a consequence of the order, in which
[`model.frame`](https://rdrr.io/r/stats/model.frame.html) applies
operations, `subset` and `na.action` would be applied after evaluation
of `r()`. That would drop some rank values from the final model frame
and returned coefficients and standard errors could no longer be
correct. The user must handle NA values and filter the data on their own
prior to usage in `lmranks`.

Wrapping `r()` with other functions (like `log(r(x))`) will not
recognize correctly the mark (because it will not be caught in
`terms(formula, specials = "r")`). The ranks will be calculated
correctly, but their transformation will be treated later in `lm` as a
regular regressor. This means that the corresponding regression
coefficient will be calculated correctly, but the standard errors,
statistics etc. will not.

`r`, `.r_predict` and `.r_cache` are special expressions, used
internally to interpret `r` mark correctly. Do not use them in
`formula`.

A number of methods defined for `lm` do not yield theoretically correct
results when applied to `lmranks` objects; errors or warnings are raised
in those instances. Also, the `df.residual` component is set to NA,
since the notion of effects of freedom for the rank models is not
theoretically established (at time of 1.2 release).

## References

Chetverikov and Wilhelm (2023), "Inference for Rank-Rank Regressions".
[arXiv preprint arXiv:2310.15512](https://arxiv.org/pdf/2310.15512)

## See also

[`lm`](https://rdrr.io/r/stats/lm.html) for details about other
arguments; [`frank`](irank.md).

Generic functions [`coef`](https://rdrr.io/r/stats/coef.html),
[`effects`](https://rdrr.io/r/stats/effects.html),
[`residuals`](https://rdrr.io/r/stats/residuals.html),
[`fitted`](https://rdrr.io/r/stats/fitted.values.html),
[`model.frame`](https://rdrr.io/r/stats/model.frame.html),
[`model.matrix`](https://rdrr.io/r/stats/model.matrix.html),
[`update`](https://rdrr.io/r/stats/update.html) .

## Examples

``` r
# rank-rank regression:
X <- rnorm(500)
Y <- X + rnorm(500)
rrfit <- lmranks(r(Y) ~ r(X))
summary(rrfit)
#> 
#> Call:
#> lmranks(formula = r(Y) ~ r(X))
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -0.65978 -0.15419 -0.00846  0.15715  0.57363 
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.15930    0.01258   12.66   <2e-16 ***
#> r(X)         0.68203    0.02511   27.16   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 

# naive version of the rank-rank regression:
RY <- frank(Y, increasing = TRUE, omega = 1)
RX <- frank(X, increasing = TRUE, omega = 1)
fit <- lm(RY ~ RX)
summary(fit)
#> 
#> Call:
#> lm(formula = RY ~ RX)
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -0.65978 -0.15419 -0.00846  0.15715  0.57363 
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  0.15930    0.01895   8.407 4.45e-16 ***
#> RX           0.68203    0.03277  20.812  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual standard error: 0.2115 on 498 degrees of freedom
#> Multiple R-squared:  0.4652, Adjusted R-squared:  0.4641 
#> F-statistic: 433.1 on 1 and 498 DF,  p-value: < 2.2e-16
#> 
# the coefficient estimates are the same as in the lmranks function, but
# the standard errors, t-values, p-values are incorrect

# support of `data` argument:
data(mtcars)
lmranks(r(mpg) ~ r(hp) + ., data = mtcars)
#> 
#> Call:
#> lmranks(formula = r(mpg) ~ r(hp) + ., data = mtcars)
#> 
#> Coefficients:
#> (Intercept)        r(hp)          cyl         disp           hp         drat  
#>    0.489054     0.094401    -0.047335     0.000655    -0.001717    -0.012264  
#>          wt         qsec           vs           am         gear         carb  
#>   -0.152140     0.036960    -0.038424     0.057301     0.072761    -0.008447  
#> 
# Same as above, but use the `hp` variable only through its rank
lmranks(r(mpg) ~ r(hp) + . - hp, data = mtcars)
#> 
#> Call:
#> lmranks(formula = r(mpg) ~ r(hp) + . - hp, data = mtcars)
#> 
#> Coefficients:
#> (Intercept)        r(hp)          cyl         disp         drat           wt  
#>   5.532e-01   -2.210e-01   -3.669e-02    9.615e-05   -1.096e-02   -1.141e-01  
#>        qsec           vs           am         gear         carb  
#>   2.944e-02   -5.069e-02    3.668e-02    7.065e-02   -2.891e-02  
#> 

# rank-rank regression with clusters:
G <- factor(rep(LETTERS[1:4], each = nrow(mtcars) / 4))
lmr <- lmranks(r(mpg) ~ r(hp):G, data = mtcars)
summary(lmr)
#> 
#> Call:
#> lmranks(formula = r(mpg) ~ r(hp):G, data = mtcars)
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -0.224762 -0.085317  0.004598  0.069257  0.249756 
#> 
#> Coefficients:
#>          Estimate Std. Error z value Pr(>|z|)    
#> GA        0.88517    0.07692  11.508  < 2e-16 ***
#> GB        1.04137    0.10402  10.012  < 2e-16 ***
#> GC        1.03120    0.05144  20.048  < 2e-16 ***
#> GD        1.04053    0.04822  21.577  < 2e-16 ***
#> r(hp):GA -0.71331    0.13239  -5.388 7.13e-08 ***
#> r(hp):GB -1.03877    0.19561  -5.311 1.09e-07 ***
#> r(hp):GC -1.07883    0.12609  -8.556  < 2e-16 ***
#> r(hp):GD -0.75260    0.08689  -8.662  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
model.matrix(lmr)
#>                     GA GB GC GD r(hp):GA r(hp):GB r(hp):GC r(hp):GD
#> Mazda RX4            1  0  0  0  0.43750  0.00000  0.00000  0.00000
#> Mazda RX4 Wag        1  0  0  0  0.43750  0.00000  0.00000  0.00000
#> Datsun 710           1  0  0  0  0.21875  0.00000  0.00000  0.00000
#> Hornet 4 Drive       1  0  0  0  0.43750  0.00000  0.00000  0.00000
#> Hornet Sportabout    1  0  0  0  0.68750  0.00000  0.00000  0.00000
#> Valiant              1  0  0  0  0.31250  0.00000  0.00000  0.00000
#> Duster 360           1  0  0  0  0.93750  0.00000  0.00000  0.00000
#> Merc 240D            1  0  0  0  0.06250  0.00000  0.00000  0.00000
#> Merc 230             0  1  0  0  0.00000  0.25000  0.00000  0.00000
#> Merc 280             0  1  0  0  0.00000  0.53125  0.00000  0.00000
#> Merc 280C            0  1  0  0  0.00000  0.53125  0.00000  0.00000
#> Merc 450SE           0  1  0  0  0.00000  0.78125  0.00000  0.00000
#> Merc 450SL           0  1  0  0  0.00000  0.78125  0.00000  0.00000
#> Merc 450SLC          0  1  0  0  0.00000  0.78125  0.00000  0.00000
#> Cadillac Fleetwood   0  1  0  0  0.00000  0.81250  0.00000  0.00000
#> Lincoln Continental  0  1  0  0  0.00000  0.84375  0.00000  0.00000
#> Chrysler Imperial    0  0  1  0  0.00000  0.00000  0.87500  0.00000
#> Fiat 128             0  0  1  0  0.00000  0.00000  0.15625  0.00000
#> Honda Civic          0  0  1  0  0.00000  0.00000  0.03125  0.00000
#> Toyota Corolla       0  0  1  0  0.00000  0.00000  0.09375  0.00000
#> Toyota Corona        0  0  1  0  0.00000  0.00000  0.28125  0.00000
#> Dodge Challenger     0  0  1  0  0.00000  0.00000  0.59375  0.00000
#> AMC Javelin          0  0  1  0  0.00000  0.00000  0.59375  0.00000
#> Camaro Z28           0  0  1  0  0.00000  0.00000  0.93750  0.00000
#> Pontiac Firebird     0  0  0  1  0.00000  0.00000  0.00000  0.68750
#> Fiat X1-9            0  0  0  1  0.00000  0.00000  0.00000  0.15625
#> Porsche 914-2        0  0  0  1  0.00000  0.00000  0.00000  0.18750
#> Lotus Europa         0  0  0  1  0.00000  0.00000  0.00000  0.46875
#> Ford Pantera L       0  0  0  1  0.00000  0.00000  0.00000  0.96875
#> Ferrari Dino         0  0  0  1  0.00000  0.00000  0.00000  0.68750
#> Maserati Bora        0  0  0  1  0.00000  0.00000  0.00000  1.00000
#> Volvo 142E           0  0  0  1  0.00000  0.00000  0.00000  0.34375
#> attr(,"assign")
#> [1] 1 1 1 1 2 2 2 2
#> attr(,"contrasts")
#> attr(,"contrasts")$G
#> [1] "contr.treatment"
#> 
# Include all columns of mtcars as usual covariates:
lmranks(r(mpg) ~ (r(hp) + .):G, data = mtcars)
#> 
#> Call:
#> lmranks(formula = r(mpg) ~ (r(hp) + .):G, data = mtcars)
#> 
#> Coefficients:
#>         GA          GB          GC          GD    r(hp):GA    r(hp):GB  
#> -2.7160344   5.3943330   1.7896255   2.6161662   0.4284687  -9.5639294  
#>   r(hp):GC    r(hp):GD      GA:cyl      GB:cyl      GC:cyl      GD:cyl  
#> -3.8243359  -0.1631026  -0.0020308   0.3431973   0.1566403  -0.0745699  
#>    GA:disp     GB:disp     GC:disp     GD:disp       GA:hp       GB:hp  
#>  0.0027112  -0.0015185  -0.0076753   0.0009933  -0.0043469   0.0279069  
#>      GC:hp       GD:hp     GA:drat     GB:drat     GC:drat     GD:drat  
#>  0.0127700  -0.0014533   0.4931073          NA  -0.1430715  -0.1490258  
#>      GA:wt       GB:wt       GC:wt       GD:wt     GA:qsec     GB:qsec  
#> -0.2639947  -0.2069544   0.5371331  -0.1647164   0.1202119  -0.2349490  
#>    GC:qsec     GD:qsec       GA:vs       GB:vs       GC:vs       GD:vs  
#> -0.0866266  -0.0250372          NA          NA          NA          NA  
#>      GA:am       GB:am       GC:am       GD:am     GA:gear     GB:gear  
#>         NA          NA          NA          NA          NA          NA  
#>    GC:gear     GD:gear     GA:carb     GB:carb     GC:carb     GD:carb  
#>         NA          NA          NA          NA          NA          NA  
#> 
```
