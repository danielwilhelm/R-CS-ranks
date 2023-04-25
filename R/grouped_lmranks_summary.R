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
    I_X <- compare(global_RX, group_RX,
                   omega=group_model$omega, na.rm=FALSE)
    I_Y <- compare(global_RY, group_RY, 
                   omega=group_model$omega, na.rm=FALSE)
    group_indicator <- as.integer(attr(object, "grouping_factor")) == j
    h_sum <- sapply(1:length(coef(group_model)), function(i){
      if(is.na(coef(group_model))[i]){
        return(rep(NA, length(group_indicator)))
      }
      proj_model <- get_projection_model(group_model, i)
      h_1 <- calculate_h_1(group_model, proj_model, group_indicator)
      h_2 <- calculate_h_2(group_model, proj_model, I_X=I_X, I_Y=I_Y)
      h_3 <- calculate_h_3(group_model, proj_model, I_X=I_X)
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
  out <- numeric(length=length(group_indicator))
  out[group_indicator] <- epsilonhat * nuhat
  return(out)
}

calculate_h_2 <- function(original_model, proj_model, I_X, I_Y){
 ineq_resids <- replace_ranks_with_ineq_indicator_and_calculate_residuals(
    original_model, I_X=I_X, I_Y=I_Y)
 # ineq_resids is N_g x N matrix
 nuhat <- resid(proj_model)
 return(colSums(ineq_resids * nuhat) / nrow(I_X))
}

calculate_h_3 <- function(original_model, proj_model, I_X){
  epsilonhat <- resid(original_model)
  if(proj_model$ranked_response){
    ineq_resids <- replace_ranks_with_ineq_indicator_and_calculate_residuals(
      proj_model, I_Y=I_X)
  } else {
    ineq_resids <- replace_ranks_with_ineq_indicator_and_calculate_residuals(
      proj_model, I_X=I_X)
  }
  return(colSums(epsilonhat * ineq_resids) / nrow(I_X))
}