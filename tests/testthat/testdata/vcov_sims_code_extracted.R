## ----setup, include=FALSE--------------------------
library(MASS)
library(AER)
library(dplyr)
library(csranks)
theta_vals <- c(0.3, 0.7)
kappa_3_vals <- 2 * sin(pi/6 * theta_vals)

# Compute corresponding k2 values to render rho = 0.5
kappa_2_vals <- 2 * sin(0.5 * asin(kappa_3_vals / 2))

# Compute k1 accordingly to acchieve different levels of wished Endogeneity
# Solve analytically when Endogeneity = 0
get_k1_zero <- function(kappa_2, kappa_3) {
  rho <- asin(kappa_2 / 2) / asin(kappa_3 / 2)
  return(2 * sin((pi / 6) * rho))
}

# Solve via optimization for other Endogeneity targets
get_k1_for_endog <- function(kappa_2, kappa_3, target_endog) {
  rho <- asin(kappa_2 / 2) / asin(kappa_3 / 2)
  
  f <- function(k1) {
    var_eps <- 1/12 + rho^2/12 - rho/pi * asin(k1/2)
    var_nu <- 1/12 + (6/pi * asin(kappa_3/2))^2/12 - (6/pi * asin(kappa_3/2))/pi * asin(kappa_3/2)
    cov_eps_nu <- 1/(2*pi) * asin(k1/2) - 1/12 * rho
    
    if (var_eps <= 0 || var_nu <= 0) return(Inf)
    
    return(abs(target_endog - cov_eps_nu / (sqrt(var_eps) * sqrt(var_nu))))
  }
  
  return(optimize(f, interval = c(-1, 1))$minimum)
}

# Construct the grid
kappa_grid <- data.frame()
for (i in 1:length(kappa_3_vals)) {
  k3 <- kappa_3_vals[i]
  k2 <- kappa_2_vals[i]
  
  k1_zero <- get_k1_zero(k2, k3)
  k1_03 <- get_k1_for_endog(k2, k3, 0.3)
  k1_07 <- get_k1_for_endog(k2, k3, 0.7)
  
  kappa_grid <- rbind(kappa_grid, 
                      data.frame(kappa_1 = c(k1_zero, k1_03, k1_07), 
                                 kappa_2 = k2, 
                                 kappa_3 = k3, 
                                 Endogeneity_Target = c(0, 0.3, 0.7)))
}

save(kappa_grid, file = "kappa_grid.RData")

knitr::opts_chunk$set(echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE)


## ----params_def------------------------------------

params <- expand.grid(n = c(500, 1000, 5000)) %>%
  merge(kappa_grid, by = NULL)  # Merge to replicate for each sample size

params$Instrument_Strength <- c(rep(0.3, 9), rep(0.7, 9))
params$mu <- list(c(2,1,0))


## ----params, echo=FALSE----------------------------
library(knitr)
param_values <- data.frame(
  Parameter = c("Sample Size (n)", "Instrument Strength", "Endogeneity"),
  Values = c("{500, 1000, 5000}", "{0.3, 0.7}", "{0, 0.3, 0.7}")
)

knitr::kable(param_values, caption = "Parameter Values Used in Simulations") 

## ----data-generation-------------------------------
gen_data <- function(n, kappa_1, kappa_2, kappa_3, mu) {
  # Define covariance matrix
  Sigma <- matrix(c(1, kappa_1, kappa_2,
                    kappa_1, 1, kappa_3,
                    kappa_2, kappa_3, 1), 
                  nrow = 3, byrow = TRUE)

  samples <- mvrnorm(n = n, mu = mu, Sigma = Sigma)
  
  # Create dataframe and transform to ranks in [0,1]
  df <- data.frame(Y = samples[,1], X = samples[,2], Z = samples[,3])
  df$rank_Y <- pnorm(df$Y - mu[1])
  df$rank_X <- pnorm(df$X - mu[2])
  df$rank_Z <- pnorm(df$Z - mu[3])
  
  df$emp_rank_Y <- rank(df$Y)/length(df$Y)
  df$emp_rank_X <- rank(df$X)/length(df$X)
  df$emp_rank_Z <- rank(df$Z)/length(df$Z)
  return(df)
}

## ----input_objects, eval=FALSE---------------------
#   object_IV <- ivreg(emp_rank_Y ~ emp_rank_X | emp_rank_Z, data = df)
#   object_SE <- lmranks(r(Y) ~ r(X), data = df)
#   object_FS <- lmranks(r(X) ~ r(Z), data = df)


## ----proj_matrix-----------------------------------
get_projection_residual_matrix <- function(object){
  regressor_dropped <- is.na(coef(object))
  if(any(regressor_dropped)){
    X <- stats::model.matrix(object)[,!regressor_dropped]
    R <- qr.R(qr(X))
  } else if(is.null(object$qr)){
    R <- qr.R(qr(stats::model.matrix(object)))
  } else {
    R <- qr.R(object$qr)
  }
  
  XTX_inv <- chol2inv(R)
  diagonal <- diag(XTX_inv)
  out <- t(t(XTX_inv) / diagonal)
  
  if(!any(regressor_dropped)){
    return(out)
  }
  full_out <- matrix(NA, nrow=length(coef(object)),
                     ncol=length(coef(object)))
  full_out[!regressor_dropped, !regressor_dropped] <- out
  full_out[regressor_dropped, !regressor_dropped] <- 0
  return(full_out)
}


## ----h1--------------------------------------------
calculate_H1 <- function(object, projection_residuals){
  original_resids <- resid(object)
  projection_residuals * original_resids
}


## ----indicator_matmult-----------------------------
ineq_indicator_matmult <- function(v, mat, omega){
  v_order <- order(v, decreasing = TRUE)
  v_ordered <- v[v_order]
  mat <- mat[v_order,,drop=FALSE]
  
  if(omega == 0){
    mat <- prepare_mat_om0(mat, v_ordered)
  } else if(omega == 1){
    mat <- prepare_mat_om1(mat, v_ordered)  
  } else {
    mat_om0 <- prepare_mat_om0(mat, v_ordered)  
    mat_om1 <- prepare_mat_om1(mat, v_ordered)
    mat <- omega*mat_om1 + (1-omega)*mat_om0
  }
  mat <- apply(mat, 2, cumsum)
  inverse_v_order <- order(v_order)
  colnames(mat) <- NULL; rownames(mat) <- NULL
  return(mat[inverse_v_order,,drop=FALSE])
}


## ----H2--------------------------------------------
calculate_H2 <- function(object_IV, object_SE, projection_residuals_FS, H1_mean=NULL){
  if(is.null(H1_mean)){
    H1_mean <- colMeans(calculate_H1(object_IV, projection_residuals_FS))
  }
  l <- get_and_separate_regressors(object_SE)
  rank_column_index <- l$rank_column_index; RX <- l$RX
  global_RX <- l$global_RX
  RY <- stats::model.response(stats::model.frame(object_SE))
  if(length(rank_column_index) > 0){
    I_X_times_proj_resids <- ineq_indicator_matmult(global_RX, projection_residuals_FS,omega=object_SE$omega)
    RX_times_proj_resids <- as.vector(global_RX %*% projection_residuals_FS)
    rowwise_rho <- get_rowwise_rho(object_IV, rank_column_index)
    delta_X_times_proj_resids <- t(I_X_times_proj_resids) - RX_times_proj_resids
    delta_X_times_proj_resids <- t(delta_X_times_proj_resids * rowwise_rho)
  } else {
    delta_X_times_proj_resids <- 0
  }
  
  if(object_SE$ranked_response){
    I_Y_times_proj_resids <- ineq_indicator_matmult(RY, projection_residuals_FS,omega=object_SE$omega)
    RY_times_proj_resids <- as.vector(RY %*% projection_residuals_FS)
    delta_Y_times_proj_resids <- t(t(I_Y_times_proj_resids) -  RY_times_proj_resids)
  }
  else
    delta_Y_times_proj_resids <- 0
  
  H2_minus_H1_mean <- (delta_Y_times_proj_resids - delta_X_times_proj_resids) / 
    stats::nobs(object_IV)
  t(t(H2_minus_H1_mean) + H1_mean)
}


## ----H3--------------------------------------------
calculate_H3 <- function(object_IV, object_SE, object_FS, projection_residual_matrix_FS, H1_mean){
  l <- get_and_separate_regressors(object_FS)
  rank_column_index <- l$rank_column_index; RZ <- l$RX
  if(length(rank_column_index)==0)
    return(0)
  global_RZ <- l$global_RX
  Z_projection_coef <- projection_residual_matrix_FS[rank_column_index,,drop=FALSE] # g columns
  original_resids <- get_original_resid_times_grouping_indicators(object_IV, object_SE)
  
  I_Z_times_orig_resids <- ineq_indicator_matmult(global_RZ, 
                                                  original_resids, 
                                                  omega=object_SE$omega) # size n x g
  RZ_times_orig_resids <- as.vector(global_RZ %*% original_resids) # size g
  delta_Z_times_orig_resids <- t(I_Z_times_orig_resids) -  RZ_times_orig_resids
  H3_minus_H1_mean <- t(delta_Z_times_orig_resids) %*% Z_projection_coef / 
    stats::nobs(object_SE)
  
  t(t(H3_minus_H1_mean) + H1_mean)
}



## ----sigma_zeta_xi, eval=FALSE---------------------
#   projection_residual_matrix_SE <- get_projection_residual_matrix(object_SE)
#   X <- stats::model.matrix(object_SE)
#   projection_residuals_SE <- X %*% projection_residual_matrix_SE
#   projection_variances <- colMeans(projection_residuals_SE * projection_residuals_FS)


## ----helper_functions, echo= FALSE-----------------
  get_and_separate_regressors <- function(model){
  if(length(model$rank_terms_indices) > 1) cli::cli_abort("Not implemented yet")
  rank_column_index <- which(model$assign %in% model$rank_terms_indices)
  if(length(rank_column_index) > 0){
    RX <- stats::model.matrix(model)[,rank_column_index]
  } else {
    RX <- integer(0)
  }
  
  if(length(rank_column_index)>1){
    global_RX <- rowSums(RX)
  } else {
    global_RX <- RX
  }
  return(list(RX=RX,
              rank_column_index=rank_column_index,
              global_RX = global_RX))
}
get_rowwise_rho <- function(object, rank_column_index){
  rho <- coef(object)[rank_column_index]
  #coef_groups <- get_coef_groups(object)
  #return(rho[coef_groups])
  return(rho)
}


get_coef_groups <- function(object){
  grouping_variable_index <- get_grouping_var_index(object)
  if(length(grouping_variable_index) == 0){
    return(rep(1, length(coef(object))))
  }
  group_levels <- levels(stats::model.frame(object)[,grouping_variable_index])
  group_indices <- sapply(1:length(coef(object)), function(i){
    coef_name <- names(coef(object))[i]
    regex <- prepare_regex_capturing_grouping_var(object, i)
    matches <- regexec(regex, coef_name, perl=TRUE)
    var_values <- regmatches(coef_name, matches)
    grouping_var_value <- var_values[[1]]["group"]
    group_idx <- which(group_levels == grouping_var_value)
    return(group_idx)
  })
  return(group_indices)
}


prepare_regex_capturing_grouping_var <- function(object, i){
  variable_table <- attr(stats::terms(object), "factors")
  grouping_variable_index <- get_grouping_var_index(object)
  var_table_column <- variable_table[,object$assign[i]]
  grouping_var_local_index <- which(var_table_column != 0) == grouping_variable_index
  
  var_names <- rownames(variable_table)[var_table_column != 0]
  var_names <- escape_special_characters(var_names)
  var_names[grouping_var_local_index] <- paste0(var_names[grouping_var_local_index], "(?<group>.*)")
  var_names[!grouping_var_local_index] <- paste0(var_names[!grouping_var_local_index], ".*")
  regex <- paste(var_names, collapse = ":")
  regex <- paste0("^", regex, "$")
  return(regex)
}
escape_special_characters <- function(v){
  regex_special_characters <- c("\\","$", "(", ")", "*", "+", ".", "?", "[", "^", "{", "|")
  escaped_v <- v
  for(spec_character in regex_special_characters){
    escaped_v <- gsub(spec_character, paste0("\\", spec_character),
                      escaped_v, fixed=TRUE)
  }
  return(escaped_v)
}

#' @noRd
get_original_resid_times_grouping_indicators <- function(object_IV, object_SE){
  original_resids <- resid(object_IV)
  #grouping_var_index <- get_grouping_var_index(object_SE)
  #if(length(grouping_var_index) > 0){
  #  grouping_var <- as.vector(stats::model.frame(object_SE)[,grouping_var_index])
  #  return(stats::model.matrix(~original_resids:grouping_var - 1))
  #} else {
    return(matrix(original_resids, ncol=1))
  #}
}
prepare_mat_om0 <- function(mat, v){
  equal_block_sizes <- diff(findIntervalIncreasing(v, TRUE))
  # 1 if the value is unique, 0 if not the last equal, k>1 if last of k equal values
  orig_mat <- mat
  mat[-1,] <- mat[-nrow(mat),] # shift 1 row, because of 0s on diagonal of I_v
  mat[1, ] <- 0 # first row of I_v is always 0
  if(all(equal_block_sizes==1)) return(mat)
  mat[which(equal_block_sizes==0)+1,] <- 0 
  om0_eq_sums <- sapply(which(equal_block_sizes>1), function(i){
    return(colSums(orig_mat[(i-equal_block_sizes[i]+1):i,,drop=FALSE]))
  })
  mat[which(equal_block_sizes>1)+1,] <- t(om0_eq_sums)
  return(mat)
}


prepare_mat_om1 <- function(mat, v){
  equal_block_sizes <- diff(c(0, findIntervalIncreasing(v, FALSE)))
  if(all(equal_block_sizes==1)) return(mat)
  om1_eq_sums <- sapply((1:nrow(mat))[equal_block_sizes>1], function(i){
    return(colSums(mat[i:(i+equal_block_sizes[i]-1),,drop=FALSE]))
  })
  mat[equal_block_sizes==0,] <- 0
  mat[equal_block_sizes>1,] <- t(om1_eq_sums)
  return(mat)
}
findIntervalIncreasing <- function(v, left.open){
  length(v) - rev(findInterval(rev(v), rev(v), left.open = !left.open))
}

## ----simulation_run, cache=TRUE--------------------
# Run simulations
j <- 3

df <- gen_data(n = params$n[j], 
               kappa_1 = params$kappa_1[j], 
               kappa_2 = params$kappa_2[j], 
               kappa_3 = params$kappa_3[j],
               mu = params$mu[[j]] )
object_IV <- ivreg(emp_rank_Y ~ emp_rank_X | emp_rank_Z, data = df)
object_SE <- lmranks(r(Y) ~ r(X), data = df)
object_FS <- lmranks(r(X) ~ r(Z), data = df)

## ----vcov------------------------------------------

projection_residual_matrix_FS <- get_projection_residual_matrix(object_FS)
Z <- stats::model.matrix(object_FS)
projection_residuals_FS <- Z %*% projection_residual_matrix_FS

H1 <- calculate_H1(object_IV, projection_residuals_FS)
H1_mean <- colMeans(H1)

H2 <- calculate_H2(object_IV, object_SE, projection_residuals_FS, H1_mean)

H3 <- calculate_H3(object_IV, object_SE, object_FS, projection_residual_matrix_FS, H1_mean)

projection_residual_matrix_SE <- get_projection_residual_matrix(object_SE)
X <- stats::model.matrix(object_SE)
projection_residuals_SE <- X %*% projection_residual_matrix_SE

projection_variances <- colMeans(projection_residuals_SE * projection_residuals_FS)
psi <- t(t(H1 + H2 + H3) / projection_variances)

sigmahat <- (t(psi) %*% psi) / (nrow(psi) ^ 2)
colnames(sigmahat) <- names(coef(object_SE))
rownames(sigmahat) <- colnames(sigmahat)
  
var_est <- sigmahat

save(df, H1, H2, H3, projection_variances, var_est, file="ivregranks_vcov_sims.rda")