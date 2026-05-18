#' @describeIn ivregranks Summary and Inference Methods for \code{"ivregranks"}
#' Objects
#'
#' @inheritParams ivreg::summary.ivreg
#' @param object An object of class \code{"ivregranks"}.
#' @param diagnostics currently not supported.
#'
#' @export
summary.ivregranks <- function(object, vcov. = NULL, df = NULL,
                               diagnostics = NULL, ...) {
  if (!is.null(vcov.)) {
    cli::cli_abort("{.var vcov.} argument is not yet supported. ")
  }
  if (!is.null(df)) {
    cli::cli_abort("{.var df} argument is not yet supported. ")
  }

  object$df.residual <- stats::nobs(object) - length(coef(object))
  outcome <- NextMethod()
  object$df.residual <- NA


  outcome$sigma <- NA
  outcome$r.squared <- NA
  outcome$adj.r.squared <- NA
  outcome$waldtest[] <- NA
  outcome$diagnostics <- NULL

  cov_matrix <- outcome$vcov
  outcome$coefficients[, 2] <- sqrt(diag(cov_matrix))
  outcome$coefficients[, 3] <- outcome$coefficients[, 1] /
    outcome$coefficients[, 2]
  outcome$coefficients[, 4] <- 2 * stats::pnorm(-abs(outcome$coefficients[, 3]))

  colnames(outcome$coefficients)[3:4] <- c("z value", "Pr(>|z|)")

  class(outcome) <- c("summary.ivregranks", class(outcome))

  return(outcome)
}

#' @rdname ivregranks
#' @export
print.summary.ivregranks <- function(x, ...) {
  x$r.squared <- x$adj.r.squared <- 0
  text <- utils::capture.output(NextMethod())
  x$r.squared <- NA
  x$adj.r.squared <- NA

  text <- text[!grepl("^Residual standard error", text)]
  text <- text[!grepl("^Multiple R-Squared", text)]
  text <- text[!grepl("^Wald test", text)]
  text <- as.list(text)
  text$sep <- "\n"
  concatenated_text <- do.call(paste, text)
  cat(concatenated_text)
  return(invisible(x))
}

#' @rdname ivregranks
#'
#' @inheritParams ivreg::confint.ivreg
#' @export
confint.ivregranks <- function(
  object, parm, level = 0.95,
  component = c("stage2", "stage1"), complete = TRUE, vcov. = NULL,
  df = NULL, ...
) {
  if (!is.null(vcov.)) {
    cli::cli_abort("{.var vcov.} argument is not yet supported. ")
  }
  if (!is.null(df)) {
    cli::cli_abort("{.var df} argument is not yet supported. ")
  }

  if (missing(parm)) {
    NextMethod(
      object = object, level = level, component = component,
      complete = complete, vcov. = vcov., ...
    )
  } else {
    NextMethod()
  }
}

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
  ## default: stage 2
  # For Z the H1 would be sth like
  # resid_Z <- resid(lm(Z ~ W)) (Bleeh)
  # resid_Y <- resid(object)
  # resid_Z * resid_Y


  regressor_dropped_fs <- is.na(coef(object))
  if (any(regressor_dropped_fs)) {
    X <- stats::model.matrix(object, component = "regressors")
    R <- qr.R(qr(X[, !regressor_dropped_fs]))
  } else {
    R <- qr.R(qr(X))
  }
  projection_residual_matrix_fs <- calculate_projection_residual_matrix(
    R,
    regressor_dropped_fs,
    length(coef(object, component = "stage1"))
  )
  projection_residuals_fs <- Z %*% projection_residual_matrix_fs

  # Why first stage? calculate_H1 would not accept this signature normally
  H1 <- calculate_H1(object, projection_residuals_fs)
  H1_mean <- colMeans(H1)
  H2 <- calculate_H2(object, projection_residuals_fs, H1_mean)
  H3 <- calculate_H3(object, projection_residual_matrix_fs, H1_mean)

  projection_residual_matrix <- get_projection_residual_matrix(object)
  X <- stats::model.matrix(object, component = "regressors")
  projection_residuals <- X %*% projection_residual_matrix
  projection_variances <- colMeans(projection_residuals *
    projection_residuals_fs)
  psi <- t(t(H1 + H2 + H3) / projection_variances)

  sigmahat <- (t(psi) %*% psi) / (nrow(psi)^2)
  colnames(sigmahat) <- names(coef(object, component = "stage2"))
  rownames(sigmahat) <- colnames(sigmahat)

  if (!complete) {
    sigmahat <- sigmahat[
      !is.na(coef(object, component = "stage2")),
      !is.na(coef(object, component = "stage2"))
    ]
  }

  return(sigmahat)
}

get_projection_residual_matrix_ivregranks <- function(object){
  stage_1_coefficients <- coef(object, component="stage1")
  regressor_dropped <- is.na(stage_1_coefficients)
  n_coef <- length(stage_1_coefficients)
  if (any(regressor_dropped)) {
    Z <- stats::model.matrix(object, component="stage1")[, !regressor_dropped]
    R <- qr.R(qr(Z))
  } else if (is.null(object[["qr1"]])) {
    R <- qr.R(qr(stats::model.matrix(object, component="stage1")))
  } else {
    R <- qr.R(object[["qr1"]])
  }
  
  calculate_projection_residual_matrix(R, regressor_dropped) 
}

calculate_projection_residual_matrix_ivregranks <- function(object){
  # Per eqn 8 from doc 'Inference for Rank-Rank Regressions with Instrumental Variables'
  # We're interested in projecting exogenous variable W_l on R(X) and other W's
  # However, on R(X) we project using 2SLS estimation.
  # We want to do that for all W's

  # For each l:
    # Step 1: project X on Z and W_-l.
    # We really don't want to do that naively (fitting a model from scratch).
    # Fortunately, we can project W_l on Z and W_-l easily using same trick as in lmranks.
  stage_1_coefficients <- coef(object, component="stage1")
  regressor_dropped <- is.na(stage_1_coefficients)

  exogenous_residual_matrix <- get_projection_residual_matrix_ivregranks(object)

  # And now the FWL theorem will convince us that we can retrieve
  # the coefficients of regression X ~ Z + W_-l by:
  # 'substituting' the linear formula for W_l from other exogeneous regressors
  # into respective coefficient.

  # unoptimized:
  # for(i in 1:nrow(exogenous_residual_matrix)){
  #   if(<i corresponds to Z>){do nothing}
  #   residual_coefficients <- exogeneous_residual_matrix[,i]

  #   W_l_coef_in_X <- stage_1_coefficients[!regressor_dropped][i]
  #   substitute <- W_l_coef_in_X * residual_coefficients * -1
  #   X_coef_without_W_l <- W_l_coef_in_X[-i] + substitute
  # }
  # opzimized:
  instrument_index <- which((1:length(stage_1_coefficients) == object[["rank_instruments_indices"]])[!regressor_dropped])
  exogeneous_residual_matrix_without_Z_projection <- exogeneous_residual_matrix[-instrument_index]

  stage_1_coefficients_cleaned <- matrix(stage_1_coefficients[!regressor_dropped][-instrument_index], nrow=1)

  substitute <- exogeneous_residual_matrix_without_Z_projection * stage_1_coefficients_cleaned * -1

  X_coefs_without_W_l <- matrix(stage_1_coefficients_cleaned, ncol=1) + substitute
  
  # suppose we have p regressors W
  # X_coefs_without_W_l is a matrix (p+1) x p
  # where in lth we have coefficients for regression
  # X ~ Z + W_-l (W_l is kept with 0)

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
  # AKA instrument_index R variable

  # Inefficient:

  c_k_div_a_k <- exogeneous_residual_matrix[instrument_index,-instrument_index] / X_coefs_without_W_l[instrument_index,]
  a_minus_e_k <- X_coefs_without_W_l - matrix(1:nrow(X_coefs_without_W_l) == instrument_index, ncol=1)
  update <- matrix(c_k_div_a_k, nrow=1) * a_minus_e_k * -1

  old_coefficients <- exogeneous_residual_matrix
  new_coefficients_l <- exogeneous_residual_matrix[-l,l] + update
}

#' Calculate H1 component for covariance estimation
#'
#' Originally defined as h_1(x, y, z) = (R_Y(y) - rhoR_X(x) - Wbeta)(R_Z(z) - Wgamma)
#'
#' @return n x p matrix
#' @noRd
calculate_H1.ivregranks <- function(object, projection_residuals, ...) {
  # We need something like
  main_model_residuals <- object$residuals
  Z_projection_residuals <- # requires a 'new model' or to be read from projection_residuals


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
