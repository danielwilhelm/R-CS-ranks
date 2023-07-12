#' Summarizing fits of linear models for ranks
#' 
#' Summary method for class "\code{lmranks}". It returns theoretically valid standard
#' errors, in comparison to naively running \code{summary(lm(...))}.
#' 
#' @param object A \code{lmranks} object.
#' @inheritParams stats::summary.lm
#' @export
summary.lmranks <- function(object, correlation = FALSE, symbolic.cor = FALSE, ...){
  if(symbolic.cor){
    cli::cli_abort("{.var symbolic.cor} are not yet implemented for {.class lmranks}.")
  }
  # call summary.lm
  object$df.residual <- nrow(object$model) - ncol(object$model)
  outcome <- NextMethod()
  object$df.residual <- NA
  
  # Mark what is unknown (for now)
  outcome$coefficients[,2:4] <- NA
  outcome$sigma <- NA
  # This one causes errors in print.summary.lm
  # If needed, we could Ctrl-C Ctrl-V and adapt the method
  #outcome$df <- c(NA, nrow(object$model) - ncol(object$model), NA)
  outcome$fstatistic <- NULL
  # Remember to handle this, once fstatistic is known
  outcome$adj.r.squared <- NA
  cov_matrix <- vcov(object, complete=FALSE)
  outcome$cov.unscaled <- matrix(NA, nrow = nrow(outcome$cov.unscaled),
                                 ncol = ncol(outcome$cov.unscaled))
  
  outcome$coefficients[, 2] <- sqrt(diag(cov_matrix))
  outcome$coefficients[, 3] <- outcome$coefficients[, 1] / outcome$coefficients[, 2]
  outcome$coefficients[, 4] <- 2*stats::pnorm(-abs(outcome$coefficients[, 3]))
  
  if(correlation)
    outcome$correlation <- stats::cov2cor(cov_matrix)
  cli::cli_warn(c("The number of residual degrees of freedom is not correct.", 
                "Also, z-value, not t-value, since the distribution used for p-value calculation is standard normal."))
  class(outcome) <- c("summary.lmranks", class(outcome))
  outcome
}

#' @export
confint.lmranks <- function(object, parm, level = 0.95, ...){
  # As is the case with confint.lm, this method returns *marginal* CIs for coefficients
  # not simultaneous
  if(missing(parm))
    stats::confint.default(object=object, level=level, ...)
  else
    stats::confint.default(object=object, parm=parm, level=level, ...)
}

#' @describeIn lmranks Calculate Variance-Covariance Matrix for a Fitted \code{lmranks} object
#' 
#' Returns the variance-covariance matrix of the regression coefficients 
#' (main parameters) of a fitted \code{lmranks} object.
#' 
#' @param complete logical indicating if the full variance-covariance matrix 
#' should be returned also in case of an over-determined system where 
#' some coefficients are undefined and \code{coef(.)} contains NAs correspondingly. 
#' When \code{complete = TRUE}, \code{vcov()} is compatible with \code{coef()} also in this singular case.
#' @importFrom stats vcov
#' @export
vcov.lmranks <- function(object, complete = TRUE, ...){
  projection_residual_matrix <- get_projection_residual_matrix(object)
  X <- stats::model.matrix(object)
  projection_residuals <- X %*% projection_residual_matrix
  original_resids <- resid(object)
  H1 <- calculate_H1(object, projection_residuals)
  
  H2 <- calculate_H2(object, projection_residuals)
  
  H1_mean <- colMeans(H1)
  H3 <- calculate_H3(object, projection_residual_matrix, H1_mean)
  
  projection_variances <- apply(projection_residuals, 2, var)
  psi <- t(t(H1 + H2 + H3) / projection_variances)
  
  sigmahat <- (t(psi) %*% psi) / (nrow(psi) ^ 2)
  colnames(sigmahat) <- names(coef(object))
  rownames(sigmahat) <- colnames(sigmahat)
  if(!complete){
    sigmahat <- sigmahat[!is.na(coef(object)),
                         !is.na(coef(object))]
  }
  return(sigmahat)
}

#' Calculate matrix giving projection residuals
#' 
#' Projections are linear models where one of X's columns is a response in terms
#' of the remaining columns.
#' 
#' @return Matrix M s.t.
#' M[i,j] = negative ith coefficient in jth projection if i != j
#' M[i,j] = 1 if i == j
#' 
#' Note, that X %*% M gives n x p matrix with residuals of jth projection in jth column.
#' 
#' Turns out, that M is closely related to V=(X^T %*% X)⁻¹:
#' M = V / diag(V), division row-wise. Proof via block matrix inverse.
#' @noRd
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

#' Calculate H1 component for covariance estimation
#' 
#' Originally defined as h_1(x,y) = (R_Y(Y)-rhoR_X(X)-Wbeta)(R_X(X) - Wgamma)
#' 
#' @return n x p matrix
#' @noRd

calculate_H1 <- function(object, projection_residuals){
  original_resids <- resid(object)
  projection_residuals * original_resids
}

#' Calculate H2 component for covariance estimation
#' 
#' Originally defined as h_2(x,y) = E[(I(y,Y)-rhoI(x,X)-Wbeta)(R_X(X) - Wgamma)]
#' Estymator in matrix notation:
#' (I_Y-rhoI_X-(Wbeta)') %*% (R_X(X)-Wgamma) / n
#' Equal to
#' I_Y %*% (R_X(X)-Wgamma) / n - 
#' rho \* I_X %*% (R_X(X)-Wgamma) / n - 
#' (Wbeta)' %*% (R_X(X)-Wgamma) / n
#' @noRd
calculate_H2 <- function(object, projection_residuals){
  l <- get_and_separate_regressors(object)
  rank_column_index <- l$rank_column_index; RX <- l$RX
  RY <- stats::model.response(stats::model.frame(object))
  if(length(rank_column_index) > 0){
    I_X_times_proj_resid <- ineq_indicator_matmult(RX, projection_residuals,omega=object$omega)
    rho <- coef(object)[rank_column_index]
    non_rank_predictor <- stats::fitted.values(object) - rho*RX
  }
  else {
    I_X_times_proj_resid <- 0
    rho <- 0
    non_rank_predictor <- stats::fitted.values(object)
  }
  
  if(object$ranked_response)
    I_Y_times_proj_resid <- ineq_indicator_matmult(RY, projection_residuals,omega=object$omega)
  else
    I_Y_times_proj_resid <- as.vector(RY %*% projection_residuals)
  
  predictor_times_proj_resid <- as.vector(non_rank_predictor %*% projection_residuals)
  out <- I_Y_times_proj_resid - rho*I_X_times_proj_resid
  out <- t(t(out) - predictor_times_proj_resid)
  return(out / stats::nobs(object))
}

#' Calculate H3 component for covariance estimation
#' 
#' Originally defined as h_3(x) = E[(R_Y(Y)-rhoR_X(X)-Wbeta)(I_X(x,X) - Wgamma)];
#' The second component depends on which projection model is considered
#' 
#' Estimator in matrix notation:
#' h_3(x) = (R_Y(Y)-rhoR_X(X)-Wbeta)' %*% [I_(x,X); W] %*% R_S / n
#' Where R_S is the projection residual matrix.
#' 
#' For a given x this higly resembles colMeans(H1).
#' The difference H3(x) - colMeans(H1)is
#'  (R_Y(Y)-rhoR_X(X)-Wbeta)'%*%(I_X(x,X) - RX)%*%R_S[r,] / n
#' (last element is a row vector from R_S matrix corresponding to ranked regressor)
#'
#' @noRd
calculate_H3 <- function(object, projection_residual_matrix, H1_mean){
  l <- get_and_separate_regressors(object)
  rank_column_index <- l$rank_column_index; RX <- l$RX
  if(length(rank_column_index)==0)
    return(0)
  
  X_projection_coef <- projection_residual_matrix[rank_column_index,]
  original_resids <- resid(object)
  I_X_times_orig_resids <- as.vector(ineq_indicator_matmult(RX, 
                                                matrix(original_resids, ncol=1),
                                                omega=object$omega)) # length n
  RX_times_orig_resids <- as.vector(RX %*% original_resids) # length 1
  delta_X_times_orig_resids <- I_X_times_orig_resids -  RX_times_orig_resids
  H3_minus_H1_mean <- delta_X_times_orig_resids %o% X_projection_coef / 
    stats::nobs(object)
  
  t(t(H3_minus_H1_mean) + H1_mean)
}

#' Extract regressors from a model object and separate rank- from usual ones
#' 
#' @return a list with entries:
#' - RX: vector of ranks of ranked regressor. May be empty.
#' - rank_column_index: which column in model.matrix corresponds to ranked regressor?
#' @noRd
get_and_separate_regressors <- function(model){
  if(length(model$rank_terms_indices) > 1) cli::cli_abort("Not implemented yet")
  rank_column_index <- which(model$assign %in% model$rank_terms_indices)
  if(length(rank_column_index) > 1) cli::cli_abort("Not implemented yet")
  if(length(rank_column_index) > 0){
    RX <- stats::model.matrix(model)[,rank_column_index]
  } else {
    RX <- integer(0)
  }
  return(list(RX=RX,
              rank_column_index=rank_column_index))
}

#' @noRd
get_global_nobs <- function(n_lequal_lesser_1, n_lequal_lesser_2){
  if(!is.null(n_lequal_lesser_1)){
    return(length(n_lequal_lesser_1$n_lequal))
  } else if(!is.null(n_lequal_lesser_2)){
    return(length(n_lequal_lesser_2$n_lequal))
  } else {
    cli::cli_abort("All regressors and responses are ordinary.")
  }
}

#' @return Numeric vector; the linear predictor part calculated from non-rank regressors.
#' @noRd
extract_nonrank_predictor <- function(model){
  l <- get_and_separate_regressors(model)
  W <- l$W; rank_column_index <- l$rank_column_index
  has_ranked_regressors <- length(rank_column_index) > 0
  if(has_ranked_regressors)
    betahat <- coef(model)[-rank_column_index]
  else
    betahat <- coef(model)
  
  # in singular fit case, some coefficients are NA
  predictor <- as.vector(W[,!is.na(betahat), drop=FALSE] %*% betahat[!is.na(betahat)])
  return(predictor)
}

#' Calculate inequality indicators and multiply with a matrix
#' 
#' @param v numeric vector
#' @param mat A matrix s.t. nrow(mat) == length(v)
#' @param omega single number
#' Inequality indicator: A matrix I_v, s.t.
#' I_v[i,j] = 1 if v[i] < v[j];
#' I_v[i,j] = omega if v[i] == v[j]; and
#' I_v[i,j] = 0 if v[i] > v[j].
#' 
#' @return I_v %*% mat
#' 
#' The trick is we do not have to do this naively (first calc I_v, then multiply).
#' And we don't want to, cause 
#' a) I_v needs O(n²) memory and 
#' b) Whole operation would take O(n²p) compute, but we repeat a lot of (simple) calculations.
#' 
#' Start with simple case: v is ordered decreasingly and has no repeating elements.
#' Then I_v is an lower triangular matrix with ones below the diagonal, 
#' omegas on diagonal and zeroes above.
#' I_v %*% mat is equivalent to taking cumsum() columnwise,
#' and then correcting for omegas on the diagonal: 
#' for C in columns of mat: 
#'    C=cumsum(omega*C + (1-omega)*c(0, C[-n]))
#' 
#' Now for a single column instead of doing O(n²) operations, we do O(n). 
#' 
#' For more general case(v ordered, with duplicates), we can use the definition of I:
#' I(a,b) = omega*i(a<=b) + (1-omega)*i(a<b)
#' And prepare the `mat` by summing entries corresponding to equal values in v.
#'
#' Finally, if v is not ordered, all we have to do is 
#' 1) permute the v and rows of mat with sorting (decrasingly) permutation of v
#' 2) proceed as in former case
#' 3) permute the rows of the result with *inverse* of sorting permutation of v
#' @noRd
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

#' @param mat Matrix s.t. nrow(mat) == length(v)
#' @param v numeric vector, sorted decreasingly
#' @return Matrix M s.t.
#' apply(M,2,cumsum) == I_v %*% mat
#' Where I_v is the inequality indicator matrix for omega = 0
#' @noRd
prepare_mat_om0 <- function(mat, v){
  
  d2 <- diff(c(0, findIntervalIncreasing(v, TRUE)))
  orig_mat <- mat
  mat[-1,] <- mat[-nrow(mat),]
  mat[d2==0,] <- 0
  if(all(d2[-1]==1)) return(mat)
  om0_eq_sums <- sapply((1:nrow(mat))[d2>1], function(i){
    return(colSums(orig_mat[(i-d2[i]):(i-1),,drop=FALSE]))
  })
  mat[d2>1,] <- t(om0_eq_sums)
  return(mat)
}

#' @param mat Matrix s.t. nrow(mat) == length(v)
#' @param v numeric vector, sorted decreasingly
#' @return Matrix M s.t.
#' apply(M,2,cumsum) == I_v %*% mat
#' Where I_v is the inequality indicator matrix for omega = 1
#' @noRd
prepare_mat_om1 <- function(mat, v){
  d1 <- diff(c(0, findIntervalIncreasing(v, FALSE)))
  if(all(d1==1)) return(mat)
  om1_eq_sums <- sapply((1:nrow(mat))[d1>1], function(i){
    return(colSums(mat[i:(i+d1[i]-1),,drop=FALSE]))
  })
  mat[d1==0,] <- 0
  mat[d1>1,] <- t(om1_eq_sums)
  return(mat)
}

#' Find Interval indices
#'
#' Given a vector of non-increasing breakpoints, find the interval containing each element;
#' If i <- findIntervalIncreasing(v), for each index j in v
#' v_{i_j} >= v_j > v_{i_{j+1}}
#' @param left.open If true, the intervals are open at left and closed at right
#' @seealso [findInterval()]
#'@noRd
findIntervalIncreasing <- function(v, left.open){
  length(v) - rev(findInterval(rev(v), rev(v), left.open = !left.open))
}

#' @importFrom stats sigma
#' @export
sigma.lmranks <- function(object, ...){
  cli::cli_abort("Not theoretically developped yet.")
}