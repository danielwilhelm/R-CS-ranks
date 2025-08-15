#' Extract regressors from a model object and separate rank- from usual ones
#'
#' @return a list with entries:
#' - RX: vector of ranks of ranked regressor. May be empty.
#' - rank_column_index: which column in model.matrix corresponds to ranked regressor?
#' @noRd
get_and_separate_regressors <- function(model_matrix, rank_column_index) {
  if (length(rank_column_index) > 0) {
    RX <- model_matrix[, rank_column_index]
  } else {
    RX <- integer(0)
  }

  if (length(rank_column_index) > 1) {
    global_RX <- rowSums(RX)
  } else {
    global_RX <- RX
  }
  return(list(
    RX = RX,
    rank_column_index = rank_column_index,
    global_RX = global_RX
  ))
}

#' @noRd
get_original_resid_times_grouping_indicators <- function(object) {
  original_resids <- resid(object)
  grouping_var_index <- get_grouping_var_index(object)
  if (length(grouping_var_index) > 0) {
    grouping_var <- as.vector(stats::model.frame(object)[, grouping_var_index])
    return(stats::model.matrix(~ original_resids:grouping_var - 1))
  } else {
    return(matrix(original_resids, ncol = 1))
  }
}

#' Each column of model.matrix(object) AKA coefficient belongs to a certain group.
#' And for each group we have 1 coefficient corresponding to ranked regressor.
#' @return a numeric vector of length same as number of coefficients
#' (number of regressor variables times number of groups) with corresponding
#' ranked coefficient (element of rho vector).
#' @noRd
get_rowwise_rho <- function(object, rank_column_index) {
  rho <- coef(object)[rank_column_index]
  coef_groups <- get_coef_groups(object)
  return(rho[coef_groups])
}

get_ranked_indices <- function(object, name) {
  if (length(object[[name]]) > 1) cli::cli_abort("Not implemented yet")
  if ("lm" %in% class(object)) {
    return(which(object$assign %in% object$rank_terms_indices))
  } else if ("ivreg" %in% class(object)) {
    return(object$ranked_instrument_indices)
  } else {
    cli::cli_abort("Object of type {.class object} not supported for this
      function")
  }
}

calculate_projection_residual_matrix <- function(R, regressor_dropped, n_coef) {
  XTX_inv <- chol2inv(R)
  diagonal <- diag(XTX_inv)
  out <- t(t(XTX_inv) / diagonal)

  if (!any(regressor_dropped)) {
    return(out)
  }
  full_out <- matrix(NA,
    nrow = n_coef,
    ncol = n_coef
  )
  full_out[!regressor_dropped, !regressor_dropped] <- out
  full_out[regressor_dropped, !regressor_dropped] <- 0
  return(full_out)
}

calculate_H1 <- function(object, ...) {
  UseMethod("calculate_H1")
}
#' Calculate H1 component for covariance estimation
#'
#' Originally defined as h_1(x,y) = (R_Y(Y)-rhoR_X(X)-Wbeta)(R_X(X) - Wgamma)
#'
#' @return n x p matrix
#' @noRd
calculate_H1.default <- function(object, projection_residuals) {
  original_resids <- resid(object)
  projection_residuals * original_resids
}

calculate_H2 <- function(object, ...) {
  UseMethod("calculate_H2")
}

#' Calculate H2 component for covariance estimation
#'
#' Originally defined as h_2(x,y) = E[(I(y,Y)-rhoI(x,X)-Wbeta)(R_X(X) - Wgamma)]
#' Estymator in matrix notation:
#' (I_Y-rhoI_X-(Wbeta)') %*% (R_X(X)-Wgamma) / n
#' Equal to
#' I_Y %*% (R_X(X)-Wgamma) / n -
#' rho \* I_X %*% (R_X(X)-Wgamma) / n -
#' (Wbeta)' %*% (R_X(X)-Wgamma) / n
#' @noRd
calculate_H2.default <- function(object, projection_residuals,
                                 H1_mean = NULL, ...) {
  l <- ..1
  RY <- ..2
  return(calculate_H2_core(
    object, l, RY, projection_residuals,
    H1_mean
  ))
}

calculate_H2_core <- function(object, l, RY, projection_residuals, H1_mean) {
  if (is.null(H1_mean)) {
    H1_mean <- colMeans(calculate_H1(object, projection_residuals))
  }
  rank_column_index <- l$rank_column_index
  RX <- l$RX
  global_RX <- l$global_RX
  if (length(rank_column_index) > 0) {
    I_X_times_proj_resids <- ineq_indicator_matmult(global_RX,
      projection_residuals,
      omega = object$omega
    )
    RX_times_proj_resids <- as.vector(global_RX %*% projection_residuals)
    rowwise_rho <- get_rowwise_rho(object, rank_column_index)
    delta_X_times_proj_resids <- t(I_X_times_proj_resids) - RX_times_proj_resids
    delta_X_times_proj_resids <- t(delta_X_times_proj_resids * rowwise_rho)
  } else {
    delta_X_times_proj_resids <- 0
  }

  if (object$ranked_response) {
    I_Y_times_proj_resids <- ineq_indicator_matmult(RY, projection_residuals, omega = object$omega)
    RY_times_proj_resids <- as.vector(RY %*% projection_residuals)
    delta_Y_times_proj_resids <- t(t(I_Y_times_proj_resids) - RY_times_proj_resids)
  } else {
    delta_Y_times_proj_resids <- 0
  }

  H2_minuse_H1_mean <- (delta_Y_times_proj_resids - delta_X_times_proj_resids) / stats::nobs(object)

  t(t(H2_minuse_H1_mean) + H1_mean)
}

calculate_H3 <- function(object, ...) {
  UseMethod("calculate_H3")
}

calculate_H3.default <- function(object, projection_residual_matrix, H1_mean, ...) {
  l <- ..1

  return(calculate_H3_core(object, l, projection_residual_matrix, H1_mean))
}

calculate_H3_core <- function(object, l, projection_residual_matrix, H1_mean) {
  rank_column_index <- l$rank_column_index
  RX <- l$RX
  if (length(rank_column_index) == 0) {
    return(0)
  }
  global_RX <- l$global_RX
  X_projection_coef <- projection_residual_matrix[rank_column_index, , drop = FALSE] # g columns
  original_resids <- get_original_resid_times_grouping_indicators(object)
  I_X_times_orig_resids <- ineq_indicator_matmult(global_RX, original_resids,
    omega = object$omega
  ) # size n x g
  RX_times_orig_resids <- as.vector(global_RX %*% original_resids) # size g
  delta_X_times_orig_resids <- t(I_X_times_orig_resids) - RX_times_orig_resids
  H3_minus_H1_mean <- t(delta_X_times_orig_resids) %*% X_projection_coef / stats::nobs(object)

  return(t(t(H3_minus_H1_mean) + H1_mean))
}
