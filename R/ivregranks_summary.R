#' Title
#'
#' @param object
#' @param vcov.
#' @param df
#' @param diagnostics
#' @param ...
#'
#' @return
#' @export
summary.ivregranks <- function(object, vcov. = NULL, df = NULL,
                               diagnostics = NULL, ...) {
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

#' Title
#'
#' @param x
#' @param ...
#'
#' @return
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

#' Title
#'
#' @param object
#' @param parm
#' @param level
#' @param component
#' @param complete
#' @param vcov.
#' @param df
#' @param ...
#'
#' @return
#' @export
confint.ivregranks <- function(
    object, parm, level = 0.95,
    component = c("stage2", "stage1"), complete = TRUE, vcov. = NULL,
    df = NULL, ...) {
  if (missing(parm)) {
    NextMethod(
      object = object, level = level, component = component,
      complete = complete, vcov. = vcov.
    )
  } else {
    NextMethod()
  }
}

#' Title
#'
#' @param object
#' @param component
#' @param complete
#' @param ...
#'
#' @return
#' @export
vcov.ivregranks <- function(object, component = c("stage2", "stage1"),
                            complete = TRUE, ...) {
  regressor_dropped_fs <- is.na(coef(object, component = "stage1"))
  Z <- stats::model.matrix(object, component = "instruments")
  if (any(regressor_dropped_fs)) {
    R <- qr.R(qr(Z[, !regressor_dropped_fs]))
  } else {
    R <- qr.R(qr(Z))
  }
  projection_residual_matrix_fs <- calculate_projection_residual_matrix(
    R,
    regressor_dropped_fs,
    length(coef(object, component = "stage1"))
  )
  projection_residuals_fs <- Z %*% projection_residual_matrix_fs

  H1 <- calculate_H1(object, projection_residuals_fs)
  H1_mean <- colMeans(H1)

  object_se <- object$object_se
  H2 <- calculate_H2(object, projection_residuals_fs, H1_mean)
  H3 <- calculate_H3(object, projection_residual_matrix_fs, H1_mean)

  projection_residual_matrix_se <- get_projection_residual_matrix(object_se)
  X <- stats::model.matrix(object_se)
  projection_residuals_se <- X %*% projection_residual_matrix_se
  projection_variances <- colMeans(projection_residuals_se *
    projection_residuals_fs)
  psi <- t(t(H1 + H2 + H3) / projection_variances)

  sigmahat <- (t(psi) %*% psi) / (nrow(psi)^2)
  colnames(sigmahat) <- names(coef(object_se))
  rownames(sigmahat) <- colnames(sigmahat)

  if (!complete) {
    sigmahat <- sigmahat[!is.na(coef(object_se)), !is.na(coef(object_se))]
  }

  return(sigmahat)
}

calculate_H1.ivregranks <- function(object, projection_residuals) {
  NextMethod()
}

calculate_H2.ivregranks <- function(object, projection_residuals, H1_mean = NULL) {
  rank_column_index <- get_ranked_indices(object$object_se, "rank_column_index")
  model_matrix_se <- stats::model.matrix(object$object_se)
  l <- get_and_separate_regressors(model_matrix_se, rank_column_index)
  RY <- stats::model.response(stats::model.frame(object$object_se))

  NextMethod(l = l, RY = RY)
}
calculate_H3.ivregranks <- function(object, projection_residual_matrix, H1_mean) {
  ranked_instrument_indices <- object$ranked_instrument_indices
  model_matrix <- stats::model.matrix(object)
  l <- get_and_separate_regressors(model_matrix, ranked_instrument_indices)

  NextMethod(l = l)
}
