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
