#' @export
#' @importFrom stats coef lm resid predict var
#' 
grouped_lmranks <- function(formula, data, grouping_factor, 
                    omega=0, na.rm=FALSE, ...){
  l <- process_lmranks_formula(formula)
  rank_terms_indices <- l$rank_terms_indices; ranked_response <- l$ranked_response
  original_call <- match.call() # for the final output
  if(length(rank_terms_indices) == 0 && !ranked_response){
    cli::cli_abort("Not implemented")
  }
  rank_env <- create_env_to_interpret_r_mark(omega, na.rm)
  environment(formula) <- rank_env
  model_frame <- lm(formula = formula, data = data, method = "model.frame")
  colnames(model_frame)[1] <- "response" # TODO: robustify
  splitted_data <- split(model_frame, f = grouping_factor, drop = TRUE)
  assign <- attr(model.matrix(formula, data=data), "assign")
  models <- lapply(splitted_data, function(df){
    model <- lm(response ~ ., data = df) # TODO: no intercept case
    model$rank_terms_indices <- rank_terms_indices
    model$assign <- assign
    model$omega <- omega
    model$ranked_response <- TRUE
    model
  })
  rank_column_index <- which(assign %in% rank_terms_indices)
  attr(models, "RX") <- model_frame[,rank_column_index] # -1 for intercept, +1 for response
  attr(models, "RY") <- model_frame$response
  attr(models, "grouping_factor") <- grouping_factor
  class(models) <- "grouped_lmranks"
  
  return(models)
}

#' @export
coef.grouped_lmranks <- function(object, ...){
  sapply(object, coef)
}