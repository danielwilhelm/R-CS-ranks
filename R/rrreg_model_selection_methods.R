#' @export
logLik.lmranks <- function(object, ...){
  cli::cli_abort("This method does not return correct results.")
}

#' @export
anova.lmranks <- function(object, ...){
  cli::cli_warn("This method might not return correct results.")
  NextMethod()
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

#' @export
plot.lmranks <- function(x,...){
  cli::cli_warn("Some plots are not correct.")
  NextMethod()
}

# AIC, BIC