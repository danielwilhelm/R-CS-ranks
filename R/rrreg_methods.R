# correct already:
# coerce, dummy.coef, family, formula, kappa, model.frame, model.matrix,
# nobs, print, qr, residuals, show, update, effects

# probably correct:
# alias, hatvalues

#' @export
logLik.lmranks <- function(object, ...){
  cli::cli_warn("This method might not return correct results.")
  NextMethod()
}

#' @export
proj.lmranks <- function(object, ...){
  cli::cli_warn("This method might not return correct results.")
  NextMethod()
}

#' @export
rstandard.lmranks <- function(object, ...){
  cli::cli_warn("This method might not return correct results.")
  NextMethod()
}

#' @export
simulate.lmranks <- function(object, nsim = 1, seed = NULL, ...){
  cli::cli_warn("This method might not return correct results.")
  NextMethod()
}

#' @export
slotsFromS3.lmranks <- function(object){
  cli::cli_warn("This method might not return correct results.")
  NextMethod()
}

#' @export
anova.lmranks <- function(object, ...){
  cli::cli_warn("This method might not return correct results.")
  NextMethod()
}

#' @export
case.names.lmranks <- function(object, ...){
  cli::cli_warn("This method might not return correct results")
  NextMethod()
}

#' @export
labels.lmranks <- function(object, ...){
  cli::cli_warn("This method might not return correct results")
  NextMethod()
}

#' @export
variable.names.lmranks <- function(object, ...){
  cli::cli_warn("This method might not return correct results")
  NextMethod()
}

#' @export
summary.lmranks <- function(object, correlation = FALSE, symbolic.cor = FALSE, ...){
  if(correlation || symbolic.cor){
    cli::cli_abort("{.var correlation} and {.var symbolic.cor} are not yet implemented for {.class lmranks}.")
  }
  # call summary.lm
  outcome <- NextMethod()
  outcome
}

#' @export
confint.lmranks <- function(object, parm, level = 0.95, ...){
  cli::cli_abort("This method returns incorrect results.",
                 "i" = "Calculation of covariances between coefficients is not yet theoretically clear.")
}

#' @export
vcov.lmranks <- function(object, ...){
  cli::cli_abort("This method returns incorrect results.",
                 "i" = "Calculation of covariances between coefficients is not yet theoretically clear.")
}

#' @export
predict.lmranks <- function(object, ...){
  cli::cli_abort("This method returns incorrect results.",
                 "i" = "Predicting a rank of a regressor is not clear.")
}

#' @export
cooks.distance.lmranks <- function(model, ...){
  cli::cli_abort("This method returns incorrect results.",
                 "i" = "Deleting and observation changes values of ranks in (in)dependent values. Standard *smart* method breaks down then.")
}

#' @export
dfbeta.lmranks <- function(model, ...){
  cli::cli_abort("This method returns incorrect results.",
                 "i" = "Deleting and observation changes values of ranks in (in)dependent values. Standard *smart* method breaks down then.")
}

#' @export
dfbetas.lmranks <- function(model, ...){
  cli::cli_abort("This method returns incorrect results.",
                 "i" = "Deleting and observation changes values of ranks in (in)dependent values. Standard *smart* method breaks down then.")
}

#' @export
influence.lmranks <- function(model, ...){
  cli::cli_abort("This method returns incorrect results.",
                 "i" = "Deleting and observation changes values of ranks in (in)dependent values. Standard *smart* method breaks down then.")
}

#' @export
rstudent.lmranks <- function(model, ...){
  cli::cli_abort("This method returns incorrect results.",
                 "i" = "Deleting and observation changes values of ranks in (in)dependent values. Standard *smart* method breaks down then.")
}

#' @export
extractAIC.lmranks <- function(model, ...){
  cli::cli_warn("NAs returned.",
                 "i" = "The notion of effective degrees of freedom for rank models is not clear.")
  c(NA, NA)
}

#' @export
add1.lmranks <- function(object, scope, ...){
  cli::cli_warn("NAs returned.",
                 "i" = "The notion of effective degrees of freedom for rank models is not clear.")
  NextMethod()
}

#' @export
drop1.lmranks <- function(object, scope, ...){
  cli::cli_abort("NAs returned.",
                 "i" = "The notion of effective degrees of freedom for rank models is not clear.")
  NextMethod()
}

#' @export
sigma.lmranks <- function(object, ...){
  cli::cli_warn("NA returned.",
                "i" = "The notion of effective degrees of freedom for rank models is not clear.")
  return(NA)
}

#' @export
plot.lmranks <- function(x,...){
  cli::cli_warn("Some plots are not correct.")
  NextMethod()
}

# AIC, BIC