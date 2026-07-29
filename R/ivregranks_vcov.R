#' @describeIn ivregranks Calculate Variance-Covariance Matrix for a Fitted
#' \code{ivregranks} object
#'
#' Returns the variance-covariance matrix of the regression coefficients
#' (main parameters) of a fitted \code{ivregranks} object. Its result is
#' theoretically valid and asymptotically consistent, in contrast to naively
#' running \code{vcov(ivreg(...))}.
#'
#' @param complete logical indicating if the full variance-covariance matrix
#' should be returned also in case of an over-determined system where
#' some coefficients are undefined and \code{coef(.)} contains NAs
#' correspondingly. When \code{complete = TRUE}, \code{vcov()} is compatible
#' with \code{coef()} also in this singular case.
#' @importFrom stats vcov
#' @export
vcov.ivregranks <- function(object, component = c("stage2", "stage1"),
                            complete = TRUE, ...) {
  component <- match.arg(component, c("stage2", "stage1"))
  if (component == "stage1") {
    return(vcov(object$object_fs, complete = complete, ...))
  }

  regressor_dropped_fs <- is.na(coef(object, component = "stage1"))

  # Per eqn 8 from doc 'Inference for Rank-Rank Regressions with Instrumental Variables'
  # We're interested in projecting exogenous variable W_l on R(X) and other W's
  # However, on R(X) we project using 2SLS estimation.
  # We want to do that for all W's

  # This one has W_l ~ RX^ + W_-l as well as RX^ ~ W
  projection_residual_matrix_stage_2 <- get_projection_residual_matrix_ivregranks(object, "stage2")

  # This one has W_l ~ Z + W_-l as well as Z ~ W
  projection_residual_matrix_stage_1 <- get_projection_residual_matrix_ivregranks(object, "stage1")

  # We don't care about projection of projection of RX^ on W, but rather Z ~ W
  instrument_index <- get_instrument_index_after_dropping_NAs(object)
  projection_residual_matrix_stage_2[, instrument_index] <- projection_residual_matrix_stage_1[, instrument_index]

  # Change for downstream logic reuse
  # Used to calculate residuals of regressions W_l ~ X^ + W_-l
  # from W and Z
  # As well as residuals of regression Z ~ W
  projection_residual_matrix_stage_2_in_terms_stage_1 <- substitute_coefs_change_base_to_stage_1(
    object,
    projection_residual_matrix_stage_2
  )

  Z <- stats::model.matrix(object, component = "instruments")[, !regressor_dropped_fs, drop = FALSE]
  projection_residuals_fs <- Z %*% projection_residual_matrix_stage_2_in_terms_stage_1

  H1 <- calculate_H1(object, projection_residuals_fs)
  H1_mean <- colMeans(H1)
  H2 <- calculate_H2(object, projection_residuals_fs, H1_mean)
  H3 <- calculate_H3(object, projection_residual_matrix_stage_2_in_terms_stage_1, H1_mean)

  X_on_W_coefs <- update_endogenous_coefficients_when_dropping_instrument(object, projection_residual_matrix_stage_1)
  X <- stats::model.matrix(object, "regressors")[, object$endogenous]
  X_on_W_residuals <- X - Z[, -instrument_index, drop = FALSE] %*% X_on_W_coefs
  projection_variances <- get_projection_variances(object, X_on_W_residuals, projection_residuals_fs)

  psi <- t(t(H1 + H2 + H3) / projection_variances)

  raw_sigmahat <- (t(psi) %*% psi) / (nrow(psi)^2)

  sigmahat <- postprocess_sigmahat(object, raw_sigmahat, complete)

  sigmahat
}

#' @title Get projection residual matrix for ivregranks
#' @description In calculation of covariance matrix of parameters in ivregranks,
#' we're interested in projections of exogenous variables (W) onto other exogenous variables, as well as on fitted values RX^ or instrument Z.
#' We really don't want to calculate the coefficients and residuals for these projections from scratch,
#' because that involves fitting p linear models (p - number of exogenous variables).
#' So we use the QR decomposition already calculated for purpose of second stage of main SLSS regression.
#' Turns out, the inverse of R^TR is closely related to coefficients of projection models.
#'
#' Important: this function returns a matrix with rows&columns corresponding only to
#' variables that were NOT colinear (their coefficients in second stage weren't NA).
#'
#' @return Suppose we have matrix U (exogenous variables + projected endogenous variable, in column order as in main regression).
#' Retuned is matrix P s.t. U %*% P gives residuals of projections of exogenous variables onto other exogenous variables as well as on projected endogenous variable.
#' First column corresponds to projection of first variable in U etc.
#'
#' @seealso get_projection_residual_matrix - a counterpart for lmranks.
#' @noRd
get_projection_residual_matrix_ivregranks <- function(object, component) {
  stage_to_matrix <- c("stage1" = "instruments", "stage2" = "projected")
  stage_to_qr_field <- c("stage1" = "qr1", "stage2" = "qr")
  matrix_component <- stage_to_matrix[component]
  qr_field <- stage_to_qr_field[component]
  stage_coefficients <- coef(object, component = component)
  regressor_dropped <- is.na(stage_coefficients)
  Z <- stats::model.matrix(object, component = matrix_component)[, !regressor_dropped]
  if (any(regressor_dropped) || is.null(object[[qr_field]])) {
    R <- qr.R(qr(Z))
  } else {
    R <- qr.R(object[[qr_field]])
  }

  calculate_projection_residual_matrix(R, regressor_dropped[!regressor_dropped], sum(!regressor_dropped))
}

#' Calculate coefficients of stage 1 regression with IVs
#' when dropping the regressors, each separately
#'
#' @param object an `ivregranks` object
#' @param stage_1_projection_residual_coefficients - A matrix of size pxp, where p is the
#' number of instrumental variables + number of exogenous variables.
#' ith column of the matrix relates to coefficients of regression of ith such variable in terms of everybody else,
#' s.t. the ith row is 1 and the rest are negative coefficients of respective variables.
#'
#' @return A vector of length p matrix with regression coefficients of model RX ~ W (note absence of Z).
#' @noRd
update_endogenous_coefficients_when_dropping_instrument <- function(
  object, stage_1_projection_residual_coeffiecients
) {
  # And now the FWL theorem will convince us that we can retrieve
  # the coefficients of regression X ~ W by:
  # 'substituting' the linear formula for Z from other exogeneous regressors
  # into respective coefficient.
  stage_1_coefficients <- coef(object, component = "stage1")
  regressor_dropped <- is.na(stage_1_coefficients)
  instrument_index <- get_instrument_index_after_dropping_NAs(object)

  stage_1_coefficients_cleaned <- stage_1_coefficients[!regressor_dropped]

  substitute <- stage_1_projection_residual_coeffiecients[-instrument_index, instrument_index] * stage_1_coefficients_cleaned[instrument_index] * -1

  stage_1_coefficients_cleaned[-instrument_index] + substitute
}

substitute_coefs_change_base_to_stage_1 <- function(object, projection_residual_matrix_stage_2) {
  # OK.
  # projection_residual_matrix_stage_2 has coefficients of regressions of type
  # W_l ~ W_-l + X^ (X^ comes from X~Z + W) and Z ~ W
  # and we want to express these coefficients in new base (with Z instead of X)
  # s.t. we can take model matrix from stage 1 and multiply it with outcome of this function
  # To get residuals needed for right bracket of e.g. H1 component
  instrument_index <- get_instrument_index_after_dropping_NAs(object)
  stage_1_coefs <- coef(object, component = "stage1")
  stage_1_coefs <- stage_1_coefs[!is.na(stage_1_coefs)]

  correction <- stage_1_coefs %o% projection_residual_matrix_stage_2[instrument_index, ]
  correction[, instrument_index] <- 0
  projection_residual_matrix_stage_2[instrument_index, -instrument_index] <- 0

  outcome <- projection_residual_matrix_stage_2 + correction
  outcome
}

get_projection_variances <- function(
  object,
  projection_residuals_endogenous_on_exogenous, projection_residuals_fs
) {
  instrument_index <- get_instrument_index_after_dropping_NAs(object)

  projection_residuals_ss <- projection_residuals_fs
  projection_residuals_ss[, instrument_index] <- projection_residuals_endogenous_on_exogenous

  projection_variances_in_stage_1_order <- colMeans(projection_residuals_ss * projection_residuals_fs)
  names(projection_variances_in_stage_1_order) <- names(coef(object, "stage1"))[!is.na(coef(object, "stage1"))]
  target_names <- names(coef(object, "stage2"))
  target_names[object[["endogenous"]]] <- names(coef(object, "stage1"))[instrument_index]
  target_names <- target_names[!is.na(coef(object, "stage2"))]

  projection_variances_in_stage_1_order[target_names]
}

get_instrument_index_after_dropping_NAs <- function(object) {
  # That's in model.matrix order
  stage_1_coefficients <- coef(object, component = "stage1")
  regressor_dropped <- is.na(stage_1_coefficients)

  matrix_column_corresponds_to_ranked_term <- attr(stats::model.matrix(object, "instruments"), "assign") %in% object[["rank_instruments_indices"]]
  instrument_index <- which(matrix_column_corresponds_to_ranked_term)

  instrument_index <- which((seq_along(stage_1_coefficients) == instrument_index)[!regressor_dropped])
  instrument_index
}

get_regressor_index_after_dropping_NAs <- function(object) {
  stage_2_coefficients <- coef(object, component = "stage2")
  regressor_dropped <- is.na(stage_2_coefficients)

  matrix_column_corresponds_to_ranked_term <- attr(stats::model.matrix(object, "regressors"), "assign") %in% object[["rank_terms_indices"]]
  regressor_term <- which(matrix_column_corresponds_to_ranked_term)

  regressor_index <- which((seq_along(stage_2_coefficients) == regressor_term)[!regressor_dropped])
  regressor_index
}

postprocess_sigmahat <- function(object, sigmahat, complete) {
  regressor_dropped_ss <- is.na(coef(object, component = "stage2"))
  if (complete && any(regressor_dropped_ss)) {
    full_sigmahat <- matrix(NA, nrow = length(regressor_dropped_ss), ncol = length(regressor_dropped_ss))
    full_sigmahat[!regressor_dropped_ss, !regressor_dropped_ss] <- sigmahat
    colnames(full_sigmahat) <- names(coef(object, component = "stage2"))
  } else {
    full_sigmahat <- sigmahat
    colnames(full_sigmahat) <- names(coef(object, component = "stage2"))[!regressor_dropped_ss]
  }

  rownames(full_sigmahat) <- colnames(full_sigmahat)
  full_sigmahat
}

#' Calculate H1 component for covariance estimation
#'
#' Originally defined as h_1(x, y, z) = (R_Y(y) - rhoR_X(x) - Wbeta)(R_Z(z) - Wgamma)
#'
#' @return n x p matrix
#' @noRd
#' @exportS3Method
calculate_H1.ivregranks <- function(object, projection_residuals, ...) {
  NextMethod()
}

#' Calculate H2 component for covariance estimation
#'
#' Originally defined as `h_2(x,y) = E[(I(y,Y)-rhoI(x,X)-Wbeta)(R_Z(Z) - Wgamma)]`
#' Estimator in matrix notation:
#' `(I_Y-rhoI_X-(Wbeta)') %*% (R_Z(Z)-Wgamma) / n`
#' Equal to
#' `I_Y %*% (R_Z(Z)-Wgamma) / n -`
#' `rho \* I_X %*% (R_Z(Z)-Wgamma) / n -`
#' `(Wbeta)' %*% (R_Z(Z)-Wgamma) / n`
#'
#' @noRd
#' @exportS3Method
calculate_H2.ivregranks <- function(object, projection_residuals,
                                    H1_mean = NULL, ...) {
  rank_column_index <- get_ranked_indices(
    object,
    component = "regressors"
  )
  model_matrix_seqn <- stats::model.matrix(object, component = "regressors")
  l <- get_and_separate_regressors(
    model_matrix_seqn,
    rank_column_index
  )
  RY <- stats::model.response(stats::model.frame(object))

  NextMethod(l = l, RY = RY)
}

#' Calculate H3 component for covariance estimation
#'
#' Originally defined as `h_3(x) = E[(R_Y(Y)-rhoR_X(X)-Wbeta)(I(z,Z) - Wgamma)]`;
#' The second component depends on which projection model is considered
#'
#' Estimator in matrix notation:
#' `h_3(x) = (R_Y(Y)-rhoR_X(X)-Wbeta)' %*% [I_(z,Z); W] %*% R_S / n`
#' Where R_S is the projection residual matrix.
#'
#' For a given x this higly resembles colMeans(H1).
#' The difference H3(x) - colMeans(H1)is
#' `(R_Y(Y)-rhoR_X(X)-Wbeta)'%*%(I(z,Z) - RX)%*%R_S[r,] / n`
#' (last element is a row vector from R_S matrix corresponding to ranked regressor)
#'
#' @noRd
#' @exportS3Method
calculate_H3.ivregranks <- function(object, projection_residual_matrix,
                                    H1_mean, ...) {
  rank_column_index <- get_ranked_indices(object,
    component = "instruments"
  )
  model_matrix <- stats::model.matrix(object, component = "instruments")
  l <- get_and_separate_regressors(model_matrix, rank_column_index)

  NextMethod(l = l)
}
