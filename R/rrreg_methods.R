# correct already:
# coerce, dummy.coef, family, formula, kappa, model.frame, model.matrix,
# nobs, print, qr, residuals, show, update, effects

# probably correct:
# add1, alias, drop1, hatvalues, sigma

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
  cli::cli_abort("This method returns incorrect results.")
}

vcov.lmranks <- function(object, ...){
  cli::cli_abort("This method returns incorrect results.")
}

predict.lmranks <- function(object, ...){
  cli::cli_abort("This method returns incorrect results.")
}

cooks.distance.lmranks <- function(model, ...){
  cli::cli_abort("This method returns incorrect results.")
}

dfbeta.lmranks <- function(model, ...){
  cli::cli_abort("This method returns incorrect results.")
}

dfbetas.lmranks <- function(model, ...){
  cli::cli_abort("This method returns incorrect results.")
}

influence.lmranks <- function(model, ...){
  cli::cli_abort("This method returns incorrect results.")
}

rstudent.lmranks <- function(model, ...){
  cli::cli_abort("This method returns incorrect results.")
}

extractAIC.lmranks <- function(model, ...){
  cli::cli_abort("The notion of effective degrees of freedom for rank models is not clear.")
}

plot.lmranks <- function(x,...){
  cli::cli_warn("Some plots are not correct.")
  NextMethod()
}

# AIC, BIC