#' @export
vcov.grouped_lmranks <- function(object, ...){
  lapply(object, vcov, ...)
}

#' @importFrom sandwich vcovHC
#' @export
vcovHC.grouped_lmranks <- function(x, ...){
  lapply(x, sandwich::vcovHC, ...)
}

# Some of the functions used here are defined in lmranks_summary
#' @export
calculate_grouped_lmranks_covariances <- function(object, complete = TRUE, ...){
  global_RX <- attr(object, "RX")
  global_RY <- attr(object, "RY")
  
  covariances <- lapply(1:length(object), function(j){
    group_model <- object[[j]]
    l <- get_and_separate_regressors(group_model)
    group_RX <- l$RX; rank_column_index <- l$rank_column_index
    group_RY <- model.response(model.frame(group_model))
    if(group_model$ranked_response){
      n_lequal_lesser_Y <- count_lequal_lesser(global_RY, group_RY, 
                                               return_inverse_ranking=TRUE)
    } else {
      n_lequal_lesser_Y <- NULL
    }
    if(length(group_model$rank_terms_indices) == 1){
      n_lequal_lesser_X <- count_lequal_lesser(global_RX, group_RX,
                                               return_inverse_ranking=TRUE)
    } else if(length(group_model$rank_terms_indices) == 0) {
      n_lequal_lesser_X <- NULL
    } else {
      cli::cli_abort("Not implemented")
    }
    
    group_indicator <- as.integer(attr(object, "grouping_factor")) == j
    h_sum <- sapply(1:length(coef(group_model)), function(i){
      if(is.na(coef(group_model))[i]){
        return(rep(NA, length(group_indicator)))
      }
      proj_model <- get_projection_model(group_model, i)
      h_1 <- calculate_h_1(group_model, proj_model, group_indicator)
      h_2 <- calculate_h_2(group_model, proj_model, n_lequal_lesser_X=n_lequal_lesser_X, 
                           n_lequal_lesser_Y=n_lequal_lesser_Y)
      h_3 <- calculate_h_3(group_model, proj_model, n_lequal_lesser_X=n_lequal_lesser_X)
      var_estimate <- sum(resid(proj_model)^2) / length(global_RX)
      (h_1 + h_2 + h_3) / var_estimate
    })
    sigmahat <- (t(h_sum) %*% h_sum) / (length(global_RX) ^ 2)
    # first division by n caused by E[h^2] estimate, second by CLT correction
    colnames(sigmahat) <- names(coef(group_model))
    rownames(sigmahat) <- colnames(sigmahat)
    sigmahat
  })
  
  names(covariances) <- names(object)
  return(covariances)
}

calculate_h_1 <- function(original_model, proj_model, group_indicator){
  epsilonhat <- resid(original_model)
  nuhat <- resid(proj_model)
  out <- numeric(length=length(group_indicator)) # 0 by default
  out[group_indicator] <- epsilonhat * nuhat
  return(out)
}

calculate_h_2 <- function(original_model, proj_model, n_lequal_lesser_X, n_lequal_lesser_Y){
 nuhat <- resid(proj_model)
 # Correct for mean in next function. Expectation over whole dataset, not just cluster
 nuhat <- nuhat * stats::nobs(original_model) / get_global_nobs(n_lequal_lesser_X, NULL) 
 return(calculate_weighted_ineq_resid_means(
   original_model, weights=nuhat, n_lequal_lesser_X=n_lequal_lesser_X, n_lequal_lesser_Y=n_lequal_lesser_Y
 ))
}

calculate_h_3 <- function(original_model, proj_model, n_lequal_lesser_X){
  if(is.null(n_lequal_lesser_X)){
    return(0)
  }
  epsilonhat <- resid(original_model)
  # Correct for mean in next function. Expectation over whole dataset, not just cluster
  epsilonhat <- epsilonhat * stats::nobs(original_model) / get_global_nobs(n_lequal_lesser_X, NULL) 
  if(proj_model$ranked_response){
    return(calculate_weighted_ineq_resid_means(
      proj_model, weights=epsilonhat, n_lequal_lesser_Y=n_lequal_lesser_X))
  } else {
    return(calculate_weighted_ineq_resid_means(
      proj_model, weights=epsilonhat, n_lequal_lesser_X=n_lequal_lesser_X))
  }
}