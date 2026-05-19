#' Calculate coefficients of stage 1 regression with IVs
#' when dropping the regressors, each separately
#'
#' @param object an `ivregranks` object
#' @param stage_1_projection_residual_coefficients - A matrix of size pxp, where p is the
#' number of instrumental variables + number of exogenous variables.
#' ith column of the matrix relates to coefficients of regression of ith such variable in terms of everybody else,
#' s.t. the ith row is 1 and the rest are negative coefficients of respective variables.
#'
#' @return A (pxp) matrix. Every column corresponds to regression
#' r(X) (endogeneous variable) ~ W_-l + r(Z) or r(X) ~ W (i.e. with one of regressors dropped).
#' @noRd
update_coefficients_when_dropping_regressors <- function(
  object, stage_1_projection_residual_coeffiecients
) {
  stage_1_coefficients <- coef(object, component = "stage1")
  regressor_dropped <- is.na(stage_1_coefficients)

  stage_1_coefficients_cleaned <- stage_1_coefficients[!regressor_dropped]

  substitute <- stage_1_projection_residual_coeffiecients * stage_1_coefficients_cleaned[col(stage_1_projection_residual_coeffiecients)] * -1

  stage_1_coefficients_cleaned[row(substitute)] + substitute
}
