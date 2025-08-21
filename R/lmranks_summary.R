#' @describeIn lmranks Summarizing fits of rank-rank regressions
#'
#' @param object A \code{lmranks} object.
#' @inheritParams stats::summary.lm
#' @export
summary.lmranks <- function(object, correlation = FALSE, symbolic.cor = FALSE, ...) {
  if (symbolic.cor) {
    cli::cli_abort("{.var symbolic.cor} are not yet implemented for {.class lmranks}.")
  }
  # call summary.lm
  object$df.residual <- stats::nobs(object) - length(coef(object))
  outcome <- NextMethod()
  object$df.residual <- NA

  # Mark what is unknown (for now)
  outcome$coefficients[, 2:4] <- NA
  outcome$sigma <- NA
  # This one causes errors in print.summary.lm
  # If needed, we could Ctrl-C Ctrl-V and adapt the method
  outcome$naive_df <- outcome$df
  outcome$df <- c(NA, NA, NA)
  outcome$fstatistic <- NULL
  # Remember to handle this, once fstatistic is known
  outcome$adj.r.squared <- NA
  cov_matrix <- vcov(object, complete = FALSE)
  outcome$cov.unscaled <- matrix(NA,
    nrow = nrow(outcome$cov.unscaled),
    ncol = ncol(outcome$cov.unscaled)
  )

  outcome$coefficients[, 2] <- sqrt(diag(cov_matrix))
  outcome$coefficients[, 3] <- outcome$coefficients[, 1] / outcome$coefficients[, 2]
  outcome$coefficients[, 4] <- 2 * stats::pnorm(-abs(outcome$coefficients[, 3]))

  colnames(outcome$coefficients)[3:4] <- c("z value", "Pr(>|z|)")

  if (correlation) {
    outcome$correlation <- stats::cov2cor(cov_matrix)
  }
  class(outcome) <- c("summary.lmranks", class(outcome))
  outcome
}

#' @export
print.summary.lmranks <- function(x, ...) {
  x$df <- x$naive_df
  text <- utils::capture.output(NextMethod())
  x$df <- c(NA, NA, NA)
  # remove the line about Residual standard error
  text <- text[!grepl("^Residual standard error: NA on [0-9]+ degrees of freedom$", text)]
  text <- as.list(text)
  text$sep <- "\n"
  concatenated_text <- do.call(paste, text)
  cat(concatenated_text)
  return(invisible(x))
}

#' @export
confint.lmranks <- function(object, parm, level = 0.95, ...) {
  # As is the case with confint.lm, this method returns *marginal* CIs for coefficients
  # not simultaneous
  if (missing(parm)) {
    stats::confint.default(object = object, level = level, ...)
  } else {
    stats::confint.default(object = object, parm = parm, level = level, ...)
  }
}

#' @describeIn lmranks Calculate Variance-Covariance Matrix for a Fitted \code{lmranks} object
#'
#' Returns the variance-covariance matrix of the regression coefficients
#' (main parameters) of a fitted \code{lmranks} object. Its result is theoretically valid
#' and asymptotically consistent, in contrast to naively running \code{vcov(lm(...))}.
#'
#' @param complete logical indicating if the full variance-covariance matrix
#' should be returned also in case of an over-determined system where
#' some coefficients are undefined and \code{coef(.)} contains NAs correspondingly.
#' When \code{complete = TRUE}, \code{vcov()} is compatible with \code{coef()} also in this singular case.
#' @importFrom stats vcov
#' @export
vcov.lmranks <- function(object, complete = TRUE, ...) {
  projection_residual_matrix <- get_projection_residual_matrix(object)
  X <- stats::model.matrix(object)
  projection_residuals <- X %*% projection_residual_matrix

  H1 <- calculate_H1(object, projection_residuals)
  H1_mean <- colMeans(H1)

  H2 <- calculate_H2(object, projection_residuals, H1_mean)

  H3 <- calculate_H3(object, projection_residual_matrix, H1_mean)

  projection_variances <- colMeans(projection_residuals^2)
  psi <- t(t(H1 + H2 + H3) / projection_variances)

  sigmahat <- (t(psi) %*% psi) / (nrow(psi)^2)
  colnames(sigmahat) <- names(coef(object))
  rownames(sigmahat) <- colnames(sigmahat)
  if (!complete) {
    sigmahat <- sigmahat[
      !is.na(coef(object)),
      !is.na(coef(object))
    ]
  }
  return(sigmahat)
}

#' Calculate matrix giving projection residuals
#'
#' Projections are linear models where one of X's columns is a response in terms
#' of the remaining columns.
#'
#' @return Matrix M s.t.
#' M[i,j] = negative ith coefficient in jth projection if i != j
#' M[i,j] = 1 if i == j
#'
#' Note, that X %*% M gives n x p matrix with residuals of jth projection in jth column.
#'
#' Turns out, that M is closely related to V=(X^T %*% X)⁻¹:
#' M = V / diag(V), division row-wise. Proof via block matrix inverse.
#' @noRd
get_projection_residual_matrix <- function(object) {
  regressor_dropped <- is.na(coef(object))
  n_coef <- length(coef(object))
  if (any(regressor_dropped)) {
    X <- stats::model.matrix(object)[, !regressor_dropped]
    R <- qr.R(qr(X))
  } else if (is.null(object$qr)) {
    R <- qr.R(qr(stats::model.matrix(object)))
  } else {
    R <- qr.R(object$qr)
  }

  return(calculate_projection_residual_matrix(R, regressor_dropped, n_coef))
}

#' Calculate H1 component for covariance estimation
#'
#' Originally defined as h_1(x,y) = (R_Y(Y)-rhoR_X(X)-Wbeta)(R_X(X) - Wgamma)
#'
#' @return n x p matrix
#'
#' @noRd
calculate_H1.lmranks <- function(object, projection_residuals, ...) {
  NextMethod()
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
#'
#' @noRd
calculate_H2.lmranks <- function(object, projection_residuals, H1_mean = NULL, ...) {
  rank_column_index <- get_ranked_indices(object, component = "regressors")
  model_matrix <- stats::model.matrix(object)
  l <- get_and_separate_regressors(model_matrix, rank_column_index)
  RY <- stats::model.response(stats::model.frame(object))

  NextMethod(l = l, RY = RY)
}

#' Calculate H3 component for covariance estimation
#'
#' Originally defined as h_3(x) = E[(R_Y(Y)-rhoR_X(X)-Wbeta)(I_X(x,X) - Wgamma)];
#' The second component depends on which projection model is considered
#'
#' Estimator in matrix notation:
#' h_3(x) = (R_Y(Y)-rhoR_X(X)-Wbeta)' %*% [I_(x,X); W] %*% R_S / n
#' Where R_S is the projection residual matrix.
#'
#' For a given x this higly resembles colMeans(H1).
#' The difference H3(x) - colMeans(H1)is
#'  (R_Y(Y)-rhoR_X(X)-Wbeta)'%*%(I_X(x,X) - RX)%*%R_S[r,] / n
#' (last element is a row vector from R_S matrix corresponding to ranked regressor)
#'
#' In the grouped case, the last element is a matrix with g rows, each corresponding
#' to regressor times indicator of grouping variable
#' And the left element is also a matrix with g rows, each with original residuals
#' times indicator of grouping variable
#'
#' @noRd
calculate_H3.lmranks <- function(object, projection_residual_matrix, H1_mean, ...) {
  if (length(object$rank_terms_indices) > 1) cli::cli_abort("Not implemented yet")
  rank_column_index <- which(object$assign %in% object$rank_terms_indices)
  model_matrix <- stats::model.matrix(object)
  l <- get_and_separate_regressors(model_matrix, rank_column_index)

  NextMethod(l = l)
}

#' @importFrom stats sigma
#' @export
sigma.lmranks <- function(object, ...) {
  cli::cli_abort("Not theoretically developped yet.")
}
