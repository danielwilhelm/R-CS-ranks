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

  # This one has W_l ~ Z + W_-l as well as Z ~ W
  projection_residual_matrix_stage_1 <- get_projection_residual_matrix_ivregranks(object)

  # This one has X ~ Z + W_-l, as well as X ~ W
  X_coefs_without_W_l <- update_coefficients_when_dropping_regressors(object, projection_residual_matrix_stage_1)

  # We'll use later sth like
  # X - X_coefs_without_W_l[,-instrument_index]
  # For estimation of residuals of regression X~W

  # This one has W_l ~ (X ~ W_-l + Z) + W_-l (eqn 8), as well as Z ~ W
  projection_residual_matrix_stage_2 <- calculate_projection_residual_matrix_stage_2(object, projection_residual_matrix_stage_1, X_coefs_without_W_l)

  # We'll use it for residuals of W_l ~ X + W_-l

  # Change for downstream logic reuse
  # Used to calculate residuals of regressions W_l ~ X + W_-l
  # from W_-l and Z
  # As well as residuals of regression Z ~ W
  projection_residual_matrix_stage_2_in_terms_stage_1 <- substitute_coefs_change_base_to_stage_1(
    object,
    projection_residual_matrix_stage_2
  )

  Z <- stats::model.matrix(object, component = "instruments")[, !regressor_dropped_fs]
  projection_residuals_fs <- Z %*% projection_residual_matrix_stage_2_in_terms_stage_1

  H1 <- calculate_H1(object, projection_residuals_fs)
  H1_mean <- colMeans(H1)
  H2 <- calculate_H2(object, projection_residuals_fs, H1_mean)
  H3 <- calculate_H3(object, projection_residual_matrix_stage_2_in_terms_stage_1, H1_mean)


  # Now, we need
  # - residuals from X ~ W
  # - residuals from Z ~ W
  # - residuals from W_l ~ (X~Z+W_-l) + W_-l
  projection_variances <- get_projection_variances(object, X_coefs_without_W_l, projection_residuals_fs)

  psi <- t(t(H1 + H2 + H3) / projection_variances)

  raw_sigmahat <- (t(psi) %*% psi) / (nrow(psi)^2)

  sigmahat <- postprocess_sigmahat(object, raw_sigmahat, complete)

  sigmahat
}

get_projection_residual_matrix_ivregranks <- function(object) {
  stage_1_coefficients <- coef(object, component = "stage1")
  regressor_dropped <- is.na(stage_1_coefficients)
  if (any(regressor_dropped)) {
    Z <- stats::model.matrix(object, component = "instruments")[, !regressor_dropped]
    R <- qr.R(qr(Z))
  } else if (is.null(object[["qr1"]])) {
    R <- qr.R(qr(stats::model.matrix(object, component = "instruments")))
  } else {
    R <- qr.R(object[["qr1"]])
  }

  calculate_projection_residual_matrix(R, regressor_dropped[!regressor_dropped], sum(!regressor_dropped))
}

substitute_coefs_change_base_to_stage_1 <- function(object, projection_residual_matrix_stage_2) {
  # OK.
  # projection_residual_matrix_stage_2 has coefficients of regressions of type
  # W_l ~ W_-l + X^ (X^ comes from X~Z + W_-l) and Z ~ W
  # and we want to:
  # - substitute X^ = X ~ Z + W (notet inclusion of W_l)
  # - express these coefficients in new base (with Z instead of X)
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

get_projection_variances <- function(object, X_coefs_without_W_l, projection_residuals_fs) {
  instrument_index <- get_instrument_index_after_dropping_NAs(object)
  regressor_dropped_ss <- is.na(coef(object, component = "stage2"))
  X <- model.matrix(object, component = "regressors")[, object[["endogenous"]], drop = FALSE]
  X <- X[, !regressor_dropped_ss[object[["endogenous"]]], drop = FALSE]
  W <- model.matrix(object, component = "regressors")[, object[["exogenous"]], drop = FALSE]
  W <- W[, !regressor_dropped_ss[object[["exogenous"]]], drop = FALSE]
  # That causo non-conformable arrays here:
  X_based_on_W <- W %*% X_coefs_without_W_l[-instrument_index, instrument_index, drop = FALSE]
  residuals_X_W <- X - as.vector(X_based_on_W)
  residuals_rest <- projection_residuals_fs
  projection_residuals_ss <- projection_residuals_fs
  projection_residuals_ss[, instrument_index] <- residuals_X_W

  colMeans(projection_residuals_ss * projection_residuals_fs)
}

get_instrument_index_after_dropping_NAs <- function(object) {
  stage_1_coefficients <- coef(object, component = "stage1")
  regressor_dropped <- is.na(stage_1_coefficients)

  instrument_term <- object[["instruments"]][object[["rank_instruments_indices"]]]
  instrument_index <- which((1:length(stage_1_coefficients) == instrument_term)[!regressor_dropped])
  instrument_index
}

get_regressor_index_after_dropping_NAs <- function(object) {
  stage_2_coefficients <- coef(object, component = "stage2")
  regressor_dropped <- is.na(stage_2_coefficients)

  regressor_term <- object[["rank_terms_indices"]]
  regressor_index <- which((1:length(stage_2_coefficients) == regressor_term)[!regressor_dropped])
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

calculate_projection_residual_matrix_stage_2 <- function(object, exogenous_residual_matrix, X_coefs_without_W_l) {
  # Step 2
  # So we have coefficients for:
  # X ~ Z + W_-l, call fitted values X_-l (1-dim vector), coefficients a
  # W_l ~ Z + W_-l, coefficients c
  # And we want W_l ~ X_-l + W_-l (coefficients d)
  # This is a change of vector base operation
  # Where we start from base [Z; W_-l]
  # And go to [X_-l; W_-l]
  # After some linear algebra that I leave as an exercise for the reader
  # d = c - c_k/a_k * (a - e_k)
  # where k is the index of Z in original base
  instrument_index <- get_instrument_index_after_dropping_NAs(object)

  X_coefs_without_W_l <- X_coefs_without_W_l[, -instrument_index, drop = FALSE]

  c_k_div_a_k <- exogenous_residual_matrix[instrument_index, -instrument_index] / X_coefs_without_W_l[instrument_index, ] # vector
  is_instrument_index_row <- 1:nrow(X_coefs_without_W_l) == instrument_index
  a_minus_e_k <- X_coefs_without_W_l - is_instrument_index_row[row(X_coefs_without_W_l)] # l is column-wise
  update <- c_k_div_a_k[col(a_minus_e_k)] * a_minus_e_k * -1

  old_coefficients <- exogenous_residual_matrix
  old_coefficients[, -instrument_index] <- old_coefficients[, -instrument_index] + update
  new_coefficients <- old_coefficients
  # To make it more fun, it would be perfect to return in row corresponding to Z the original coefficients (since they are exactly what we need)
  # But then we'll confuse the projections: sometimes we want to plug in estimates of RX, sometimes RZ
  # Carefully!
  new_coefficients
}

calculate_projection_residual_matrix_ivregranks <- function(object) {
  # Per eqn 8 from doc 'Inference for Rank-Rank Regressions with Instrumental Variables'
  # We're interested in projecting exogenous variable W_l on R(X) and other W's
  # However, on R(X) we project using 2SLS estimation.
  # We want to do that for all W's

  # For each l:
  # Step 1: project X on Z and W_-l.
  # We really don't want to do that naively (fitting a model from scratch).
  # Fortunately, we can project W_l on Z and W_-l easily using same trick as in lmranks.

  exogenous_residual_matrix <- get_projection_residual_matrix_ivregranks(object)

  # And now the FWL theorem will convince us that we can retrieve
  # the coefficients of regression X ~ Z + W_-l by:
  # 'substituting' the linear formula for W_l from other exogeneous regressors
  # into respective coefficient.
  X_coefs_without_W_l <- update_coefficients_when_dropping_regressors(object, exogenous_residual_matrix)

  # suppose we have p regressors W
  # X_coefs_without_W_l is a matrix (p+1) x p
  # where in lth we have coefficients for regression
  # X ~ Z + W_-l (W_l is kept with 0)

  calculate_projection_residual_matrix_stage_2(object, exogenous_residual_matrix, X_coefs_without_W_l)
}

#' Calculate H1 component for covariance estimation
#'
#' Originally defined as h_1(x, y, z) = (R_Y(y) - rhoR_X(x) - Wbeta)(R_Z(z) - Wgamma)
#'
#' @return n x p matrix
#' @noRd
calculate_H1.ivregranks <- function(object, projection_residuals, ...) {
  NextMethod()
}

#' Calculate H2 component for covariance estimation
#'
#' Originally defined as h_2(x,y) = E[(I(y,Y)-rhoI(x,X)-Wbeta)(R_Z(Z) - Wgamma)]
#' Estimator in matrix notation:
#' (I_Y-rhoI_X-(Wbeta)') %*% (R_Z(Z)-Wgamma) / n
#' Equal to
#' I_Y %*% (R_Z(Z)-Wgamma) / n -
#' rho \* I_X %*% (R_Z(Z)-Wgamma) / n -
#' (Wbeta)' %*% (R_Z(Z)-Wgamma) / n
#'
#' @noRd
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
#' Originally defined as h_3(x) = E[(R_Y(Y)-rhoR_X(X)-Wbeta)(I(z,Z) - Wgamma)];
#' The second component depends on which projection model is considered
#'
#' Estimator in matrix notation:
#' h_3(x) = (R_Y(Y)-rhoR_X(X)-Wbeta)' %*% [I_(z,Z); W] %*% R_S / n
#' Where R_S is the projection residual matrix.
#'
#' For a given x this higly resembles colMeans(H1).
#' The difference H3(x) - colMeans(H1)is
#'  (R_Y(Y)-rhoR_X(X)-Wbeta)'%*%(I(z,Z) - RX)%*%R_S[r,] / n
#' (last element is a row vector from R_S matrix corresponding to ranked regressor)
#'
#' @noRd
calculate_H3.ivregranks <- function(object, projection_residual_matrix,
                                    H1_mean, ...) {
  rank_column_index <- get_ranked_indices(object,
    component = "instruments"
  )
  model_matrix <- stats::model.matrix(object, component = "instruments")
  l <- get_and_separate_regressors(model_matrix, rank_column_index)

  NextMethod(l = l)
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
