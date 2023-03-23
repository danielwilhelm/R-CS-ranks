# correct already:
# coerce, dummy.coef, family, formula, kappa, model.frame, model.matrix,
# nobs, print, qr, residuals, show, update, effects

# probably correct:
# alias, hatvalues

logLik.lmranks <- function(object, ...){
  cli::cli_warn("This method might not return correct results.")
  NextMethod()
}

proj.lmranks <- function(object, ...){
  cli::cli_warn("This method might not return correct results.")
  NextMethod()
}

rstandard.lmranks <- function(object, ...){
  cli::cli_warn("This method might not return correct results.")
  NextMethod()
}

simulate.lmranks <- function(object, nsim = 1, seed = NULL, ...){
  cli::cli_warn("This method might not return correct results.")
  NextMethod()
}

slotsFromS3.lmranks <- function(object){
  cli::cli_warn("This method might not return correct results.")
  NextMethod()
}

anova.lmranks <- function(object, ...){
  cli::cli_warn("This method might not return correct results.")
  NextMethod()
}

case.names.lmranks <- function(object, ...){
  cli::cli_warn("This method might not return correct results")
  NextMethod()
}

labels.lmranks <- function(object, ...){
  cli::cli_warn("This method might not return correct results")
  NextMethod()
}

variable.names.lmranks <- function(object, ...){
  cli::cli_warn("This method might not return correct results")
  NextMethod()
}

summary.lmranks <- function(object, correlation = FALSE, symbolic.cor = FALSE, ...){
  if(correlation || symbolic.cor){
    cli::cli_abort("{.var correlation} and {.var symbolic.cor} are not yet implemented for {.class lmranks}.")
  }
  # call summary.lm
  outcome <- NextMethod()
  outcome
}

confint.lmranks <- function(object, parm, level = 0.95, ...){
  cli::cli_abort("This method returns incorrect results.",
                 "i" = "Calculation of covariances between coefficients is not yet theoretically clear.")
}

vcov.lmranks <- function(object, ...){
  cli::cli_abort("This method returns incorrect results.",
                 "i" = "Calculation of covariances between coefficients is not yet theoretically clear.")
}

predict.lmranks <- function(object, ...){
  cli::cli_abort("This method returns incorrect results.",
                 "i" = "Predicting a rank of a regressor is not clear.")
}

cooks.distance.lmranks <- function(model, ...){
  cli::cli_abort("This method returns incorrect results.",
                 "i" = "Deleting and observation changes values of ranks in (in)dependent values. Standard *smart* method breaks down then.")
}

dfbeta.lmranks <- function(model, ...){
  cli::cli_abort("This method returns incorrect results.",
                 "i" = "Deleting and observation changes values of ranks in (in)dependent values. Standard *smart* method breaks down then.")
}

dfbetas.lmranks <- function(model, ...){
  cli::cli_abort("This method returns incorrect results.",
                 "i" = "Deleting and observation changes values of ranks in (in)dependent values. Standard *smart* method breaks down then.")
}

influence.lmranks <- function(model, ...){
  cli::cli_abort("This method returns incorrect results.",
                 "i" = "Deleting and observation changes values of ranks in (in)dependent values. Standard *smart* method breaks down then.")
}

rstudent.lmranks <- function(model, ...){
  cli::cli_abort("This method returns incorrect results.",
                 "i" = "Deleting and observation changes values of ranks in (in)dependent values. Standard *smart* method breaks down then.")
}

extractAIC.lmranks <- function(model, ...){
  cli::cli_warn("NAs returned.",
                 "i" = "The notion of effective degrees of freedom for rank models is not clear.")
  c(NA, NA)
}

add1.lmranks <- function(object, scope, ...){
  cli::cli_warn("NAs returned.",
                 "i" = "The notion of effective degrees of freedom for rank models is not clear.")
  NextMethod()
}

drop1.lmranks <- function(object, scope, ...){
  cli::cli_abort("NAs returned.",
                 "i" = "The notion of effective degrees of freedom for rank models is not clear.")
  NextMethod()
}

sigma.lmranks <- function(object, ...){
  cli::cli_warn("NA returned.",
                "i" = "The notion of effective degrees of freedom for rank models is not clear.")
  return(NA)
}

plot.lmranks <- function(x,...){
  cli::cli_warn("Some plots are not correct.")
  NextMethod()
}

# AIC, BIC