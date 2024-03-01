#' @export
logLik.lmranks <- function(object, ...){
  cli::cli_abort("This method does not return correct results.")
}

#' @export
anova.lmranks <- function(object, ...){
  cli::cli_abort("This method is not implemented.")
}

#' @export
proj.lmranks <- function(object, ...){
  object$df.residual <- stats::nobs(object) - length(coef(object))
  out <- NextMethod()
  object$df.residual <- NA
  attr(out, "df") <- c(NA, NA, NA)
  return(out)
}

#' @export
extractAIC.lmranks <- function(fit, scale, k, ...){
  cli::cli_warn(c("NAs returned.",
                 "i" = "The notion of effective degrees of freedom for rank models is not clear."))
  c(NA, NA)
}

#' @export
add1.lmranks <- function(object, scope, ...){
  cli::cli_warn(c("NAs returned.",
                 "i" = "The notion of effective degrees of freedom for rank models is not clear."))
  NextMethod()
}

#' @export
drop1.lmranks <- function(object, scope, ...){
  cli::cli_warn(c("NAs returned.",
                 "i" = "The notion of effective degrees of freedom for rank models is not clear."))
  NextMethod()
}

# AIC, BIC