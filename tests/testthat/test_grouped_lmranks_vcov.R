###############
### lmranks ###
###############

compare_for_tests <- function(i,j,v,omega=1){
  omega * (v[i] <= v[j]) + (1-omega) * (v[i] < v[j])
}

data(mtcars)
mtcars2 <- mtcars
G <- factor(rep(LETTERS[1:4], each=nrow(mtcars) / 4))
mtcars2$G <- G

test_that("summary does not raise errors", {
  mod <- lmranks(r(mpg) ~ (r(cyl) + disp):G, data=mtcars2)
  expect_warning(summary(mod), "degrees of freedom")
})

test_that("vcov passes shallow checks", {
  model <- lmranks(r(mpg) ~ (r(cyl) + disp):G, data=mtcars2)
  V <- vcov(model)
  
  expect_true(isSymmetric(V))
  vals <- eigen(V, only.values = TRUE)
  expect_true(all(vals$values> 0))
  expect_equal(colnames(V),
               c("(Intercept)", "r(cyl)", "disp"))
  expect_equal(rownames(V),
               c("(Intercept)", "r(cyl)", "disp"))
})

test_that("vcov works for singular model matrix", {
  W <- cbind(mtcars$disp, mtcars$disp)
  mod <- lmranks(r(mpg) ~ (r(cyl) + W):G, data=mtcars2)
  cov_mat <- vcov(mod, complete=TRUE)
  
  expect_true(all(!is.na(cov_mat[1:12, 1:12])))
  expect_true(all(is.na(cov_mat[13:16,])))
  expect_true(all(is.na(cov_mat[,13:16])))
})

test_that("vcov works for singular model matrix, complete=FALSE", {
  data(mtcars)
  W <- cbind(mtcars$disp, mtcars$disp)
  mod <- lmranks(r(mpg) ~ (r(cyl) + W):G, data=mtcars2)
  cov_mat <- vcov(mod, complete=FALSE)
  
  expect_true(all(!is.na(cov_mat)))
  expect_equal(c(nrow(cov_mat), ncol(cov_mat)),
               c(12,12))
})

test_that("get_projection_residual_matrix works in grouped case", {
  X <- as.matrix(mtcars) 
  X <- X[,-((ncol(X) - 3):ncol(X))] # nonidentifiability
  XX <- model.matrix(~X:G-1)
  n_groups <- length(levels(mtcars2$G))
  expected <- matrix(0,nrow = ncol(XX),
                     ncol=ncol(XX))
  expected_resids <- matrix(0, nrow=nrow(X), ncol=ncol(XX))
  for(j in 1:ncol(X)){
    jj <- j
    for(k in 1:n_groups){
      g <- levels(mtcars2$G)[k]
      X_minusj <- X[mtcars2$G == g,-j]
      Y <- X[mtcars2$G == g,j]
      m <- lm(Y~X_minusj-1)
      coef_vector <- numeric(ncol(X))
      coef_vector[-j] <- -coef(m)
      coef_vector[j] <- 1
      group_indices <- 1:ncol(X) + ncol(X)*(k-1)
      expected[group_indices, jj] <- coef_vector
      expected_resids[mtcars2$G == g, jj] <- resid(m)
      jj <- jj + ncol(X)
    }
  }
  object <- list(qr=qr(XX))
  actual <- get_projection_residual_matrix(object)
  expect_equivalent(expected, actual)
  expect_equivalent(expected_resids, XX %*% actual)
  
  # got: p1:g1, p2:g1,...,pn:g1,p1:g2:,...
  # want: p1:g1, p1:g2,...,p1:gm,p2:g1,...
  rearranging_perm <- as.vector(matrix(1:ncol(XX), ncol=ncol(X), byrow=TRUE))
  expected_2 <- expected[rearranging_perm, rearranging_perm]
  object_2 <- list(qr=qr(XX[,rearranging_perm]))
  actual_2 <- get_projection_residual_matrix(object_2)
  expect_equivalent(expected_2, actual_2)
})

test_that("get_projection_residual_matrix works for singular matrix in grouped case", {
  data(mtcars)
  X <- as.matrix(mtcars) 
  X <- X[,-((ncol(X) - 3):ncol(X))] # nonidentifiability
  X[G=="A",2] <- X[G=="A",1]
  XX <- model.matrix(~X:G-1)
  n_groups <- length(levels(mtcars2$G))
  expected <- matrix(0,nrow = ncol(XX),
                     ncol=ncol(XX))
  expected[,2] <- NA
  for(j in 1:ncol(X)){
    jj <- j
    for(k in 1:n_groups){
      if(j==2 && k==1){
        jj <- jj + ncol(X)
        next
      }
      g <- levels(mtcars2$G)[k]
      if(k==1){
        X_minusj <- X[mtcars2$G == g,-c(2,j)]
      } else {
        X_minusj <- X[mtcars2$G == g,-j]
      }
      Y <- X[mtcars2$G == g,j]
      m <- lm(Y~X_minusj-1)
      coef_vector <- numeric(ncol(X))
      if(k==1){
        coef_vector[-c(2,j)] <- -coef(m)
        coef_vector[j] <- 1
        coef_vector[2] <- 0
      } else {
        coef_vector[-j] <- -coef(m)
        coef_vector[j] <- 1 
      }
      group_indices <- 1:ncol(X) + ncol(X)*(k-1)
      expected[group_indices, jj] <- coef_vector
      jj <- jj + ncol(X)
    }
  }
  
  Y <- mtcars[,1]
  object <- lm(Y~XX-1)
  actual <- get_projection_residual_matrix(object)
  expect_equivalent(expected, actual)
})

test_that("get_and_separate_regressors works with grouping",{
  data(mtcars)
  mtcars2 <- mtcars
  mtcars2$cyl <- factor(mtcars2$cyl)
  n_groups <- length(levels(mtcars2$cyl))
  model_1 <- lm(mpg ~ disp:cyl, data=mtcars2)
  model_1$rank_terms_indices <- 1
  expected_RX <- model.matrix(model_1)[,1:n_groups + 1] # Intercept
  expected_out_1 <- list(RX = expected_RX,
                         rank_column_index = 1:n_groups+1, # Intercept
                         global_RX = rowSums(expected_RX))
  expect_equal(get_and_separate_regressors(model_1),
               expected_out_1)
})

test_that("get_coef_groups works",{
  nasty_df <- mtcars2
  nasty_df$G <- factor(c("AB", "AC", "AD", "AE")[as.integer(G)])
  nasty_df$GA <- sample(G)
  
  m1 <- lmranks(r(mpg) ~ r(disp):G + GA:G - 1, data=nasty_df)
  m2 <- lmranks(r(mpg) ~ GA:G + r(disp):G - 1, data=nasty_df)
  m3 <- lmranks(r(mpg) ~ (r(disp) + .):G, data=nasty_df[,c("mpg", "disp", "GA", "hp", "G")])
  
  expected_1 <- c(1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4)
  expected_2 <- c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,1,2,3,4)
  expected_3 <- rep(1:4, times=7)
  
  actual <- get_coef_groups(m1)
  expect_equivalent(actual, expected_1)
  actual <- get_coef_groups(m2)
  expect_equivalent(actual, expected_2)
  actual <- get_coef_groups(m3)
  expect_equivalent(actual, expected_3)
  
})

#########################
### High-level checks ###
#########################

test_that("h1 works for ranked regressor with no covariates with grouping", {
  load(test_path("testdata", "grouped_lmranks_cov_sigmahat_covariates_FALSE.rda"))
  G <- factor(G)
  res <- lmranks(r(Y) ~ r(X):G, omega=1)
  proj_residuals <- model.matrix(res) %*% get_projection_residual_matrix(res)
  
  h1_lmranks <- calculate_H1(res, proj_residuals)
  h11_lmranks <- h1_lmranks[,3]
  h12_lmranks <- h1_lmranks[,4]
  expect_equivalent(h11_lmranks, H11)
  expect_equivalent(h12_lmranks, H12)
})

test_that("h2 works for ranked regressor with no covariates with grouping", {
  load(test_path("testdata", "grouped_lmranks_cov_sigmahat_covariates_FALSE.rda"))
  G <- factor(G)
  res <- lmranks(r(Y) ~ r(X):G, omega=1)
  proj_residuals <- model.matrix(res) %*% get_projection_residual_matrix(res)
  
  h2_lmranks <- calculate_H2(res, proj_residuals)
  h21_lmranks <- h2_lmranks[,3]
  h22_lmranks <- h2_lmranks[,4]
  expect_equivalent(h21_lmranks, H21)
  expect_equivalent(h22_lmranks, H22)
})

test_that("h3 works for ranked regressor with no covariates with grouping", {
  load(test_path("testdata", "grouped_lmranks_cov_sigmahat_covariates_FALSE.rda"))
  res <- grouped_lmranks(r(Y) ~ r(X), data=data.frame(Y=Y,X=X), grouping_factor=G, omega=1)
  
  n_lequal_lesser_X_1 <- count_lequal_lesser(X, X[G==1], 
                                             return_inverse_ranking=TRUE)
  proj_model_1 <- get_projection_model(res[[1]], 2)
  h31_lmranks <- calculate_h_3(res[[1]], proj_model_1, 
                               n_lequal_lesser_X=n_lequal_lesser_X_1)
  expect_equivalent(h31_lmranks, H31)
  
  n_lequal_lesser_X_2 <- count_lequal_lesser(X, X[G==2], 
                                             return_inverse_ranking=TRUE)
  proj_model_2 <- get_projection_model(res[[2]], 2)
  h32_lmranks <- calculate_h_3(res[[2]], proj_model_2, 
                               n_lequal_lesser_X=n_lequal_lesser_X_2)
  expect_equivalent(h32_lmranks, H32)
})

test_that("vcov produces correct asymptotic variance estimate of rank-rank slope with no covariates", {
  load(test_path("testdata", "grouped_lmranks_cov_sigmahat_covariates_FALSE.rda"))
  res <- grouped_lmranks(r(Y) ~ r(X), data=data.frame(Y=Y,X=X), grouping_factor=G, omega=1)
  cov_mat <- vcov(res)
  sigma2hat.grouped_lmranks <- c(cov_mat[2,2]*n, cov_mat[4,4]*n)
  expect_equal(sigma2hat, sigma2hat.grouped_lmranks, tolerance=1e-5)
})

test_that("h1 works for ranked regressor with covariates and grouping", {
  load(test_path("testdata", "grouped_lmranks_cov_sigmahat_covariates_TRUE.rda"))
  res <- grouped_lmranks(r(Y) ~ r(X)+W-1, data=data.frame(Y=Y,X=X), grouping_factor=G, omega=1)
  proj_model_1 <- get_projection_model(res[[1]], 1)
  g1 <- as.integer(attr(res, "grouping_factor")) == 1
  h11_lmranks <- calculate_h_1(res[[1]], proj_model_1, g1)
  expect_equivalent(h11_lmranks, H11)
  
  proj_model_2 <- get_projection_model(res[[2]], 1)
  g2 <- as.integer(attr(res, "grouping_factor")) == 2
  h12_lmranks <- calculate_h_1(res[[2]], proj_model_2, g2)
  expect_equivalent(h12_lmranks, H12)
})

test_that("h2 works for ranked regressor with covariates and grouping", {
  load(test_path("testdata", "grouped_lmranks_cov_sigmahat_covariates_TRUE.rda"))
  res <- grouped_lmranks(r(Y) ~ r(X)+W-1, data=data.frame(Y=Y,X=X), grouping_factor=G, omega=1)
  
  n_lequal_lesser_X_1 <- count_lequal_lesser(X, X[G==1], 
                                             return_inverse_ranking=TRUE)
  n_lequal_lesser_Y_1 <- count_lequal_lesser(Y, Y[G==1], 
                                             return_inverse_ranking=TRUE)
  proj_model_1 <- get_projection_model(res[[1]], 1)
  h21_lmranks <- calculate_h_2(res[[1]], proj_model_1, 
                               n_lequal_lesser_X=n_lequal_lesser_X_1, 
                               n_lequal_lesser_Y=n_lequal_lesser_Y_1)
  expect_equivalent(h21_lmranks, H21)
  
  n_lequal_lesser_X_2 <- count_lequal_lesser(X, X[G==2], 
                                             return_inverse_ranking=TRUE)
  n_lequal_lesser_Y_2 <- count_lequal_lesser(Y, Y[G==2], 
                                             return_inverse_ranking=TRUE)
  proj_model_2 <- get_projection_model(res[[2]], 1)
  h22_lmranks <- calculate_h_2(res[[2]], proj_model_2, 
                               n_lequal_lesser_X=n_lequal_lesser_X_2, 
                               n_lequal_lesser_Y=n_lequal_lesser_Y_2)
  expect_equivalent(h22_lmranks, H22)
})

test_that("h3 works for ranked regressor with covariates and grouping", {
  load(test_path("testdata", "grouped_lmranks_cov_sigmahat_covariates_TRUE.rda"))
  res <- grouped_lmranks(r(Y) ~ r(X)+W-1, data=data.frame(Y=Y,X=X), grouping_factor=G, omega=1)
  
  n_lequal_lesser_X_1 <- count_lequal_lesser(X, X[G==1], 
                                             return_inverse_ranking=TRUE)
  proj_model_1 <- get_projection_model(res[[1]], 1)
  h31_lmranks <- calculate_h_3(res[[1]], proj_model_1, 
                               n_lequal_lesser_X=n_lequal_lesser_X_1)
  expect_equivalent(h31_lmranks, H31)
  
  n_lequal_lesser_X_2 <- count_lequal_lesser(X, X[G==2], 
                                             return_inverse_ranking=TRUE)
  proj_model_2 <- get_projection_model(res[[2]], 1)
  h32_lmranks <- calculate_h_3(res[[2]], proj_model_2, 
                               n_lequal_lesser_X=n_lequal_lesser_X_2)
  expect_equivalent(h32_lmranks, H32)
})


test_that("vcov produces correct asymptotic variance estimate of rank-rank slope with covariates", {
  load(test_path("testdata", "grouped_lmranks_cov_sigmahat_covariates_TRUE.rda"))
  res <- grouped_lmranks(r(Y) ~ r(X) + W - 1, data=data.frame(Y=Y,X=X), grouping_factor=G, omega=1)
  cov_mat <- vcov(res)
  sigma2hat.grouped_lmranks <- c(cov_mat[1,1]*n, cov_mat[5,5]*n)
  expect_equal(sigma2hat, sigma2hat.grouped_lmranks, tolerance=1e-5)
  
  W_no_intercept <- W[,-1]
  res <- grouped_lmranks(r(Y) ~ r(X) + W_no_intercept, data=data.frame(Y=Y,X=X), grouping_factor=G, omega=1)
  cov_mat <- vcov(res)
  sigma2hat.grouped_lmranks <- c(cov_mat[2,2]*n, cov_mat[6,6]*n)
  expect_equal(sigma2hat, sigma2hat.grouped_lmranks, tolerance=1e-5)
})

#######################
### grouped_lmranks ###
#######################

test_that("h1 works for ranked regressor with no covariates with grouping", {
  load(test_path("testdata", "grouped_lmranks_cov_sigmahat_covariates_FALSE.rda"))
  res <- grouped_lmranks(r(Y) ~ r(X), data=data.frame(Y=Y,X=X), grouping_factor=G, omega=1)
  proj_model_1 <- get_projection_model(res[[1]], 2)
  g1 <- as.integer(attr(res, "grouping_factor")) == 1
  h11_lmranks <- calculate_h_1(res[[1]], proj_model_1, g1)
  expect_equivalent(h11_lmranks, H11)
  
  proj_model_2 <- get_projection_model(res[[2]], 2)
  g2 <- as.integer(attr(res, "grouping_factor")) == 2
  h12_lmranks <- calculate_h_1(res[[2]], proj_model_2, g2)
  expect_equivalent(h12_lmranks, H12)
})

test_that("h2 works for ranked regressor with no covariates with grouping", {
  load(test_path("testdata", "grouped_lmranks_cov_sigmahat_covariates_FALSE.rda"))
  res <- grouped_lmranks(r(Y) ~ r(X), data=data.frame(Y=Y,X=X), grouping_factor=G, omega=1)
  
  n_lequal_lesser_X_1 <- count_lequal_lesser(X, X[G==1], 
                                           return_inverse_ranking=TRUE)
  n_lequal_lesser_Y_1 <- count_lequal_lesser(Y, Y[G==1], 
                                             return_inverse_ranking=TRUE)
  proj_model_1 <- get_projection_model(res[[1]], 2)
  h21_lmranks <- calculate_h_2(res[[1]], proj_model_1, 
                               n_lequal_lesser_X=n_lequal_lesser_X_1, 
                               n_lequal_lesser_Y=n_lequal_lesser_Y_1)
  expect_equivalent(h21_lmranks, H21)
  
  n_lequal_lesser_X_2 <- count_lequal_lesser(X, X[G==2], 
                                             return_inverse_ranking=TRUE)
  n_lequal_lesser_Y_2 <- count_lequal_lesser(Y, Y[G==2], 
                                             return_inverse_ranking=TRUE)
  proj_model_2 <- get_projection_model(res[[2]], 2)
  h22_lmranks <- calculate_h_2(res[[2]], proj_model_2, 
                               n_lequal_lesser_X=n_lequal_lesser_X_2, 
                               n_lequal_lesser_Y=n_lequal_lesser_Y_2)
  expect_equivalent(h22_lmranks, H22)
})

test_that("h3 works for ranked regressor with no covariates with grouping", {
  load(test_path("testdata", "grouped_lmranks_cov_sigmahat_covariates_FALSE.rda"))
  res <- grouped_lmranks(r(Y) ~ r(X), data=data.frame(Y=Y,X=X), grouping_factor=G, omega=1)
  
  n_lequal_lesser_X_1 <- count_lequal_lesser(X, X[G==1], 
                                             return_inverse_ranking=TRUE)
  proj_model_1 <- get_projection_model(res[[1]], 2)
  h31_lmranks <- calculate_h_3(res[[1]], proj_model_1, 
                               n_lequal_lesser_X=n_lequal_lesser_X_1)
  expect_equivalent(h31_lmranks, H31)
  
  n_lequal_lesser_X_2 <- count_lequal_lesser(X, X[G==2], 
                                             return_inverse_ranking=TRUE)
  proj_model_2 <- get_projection_model(res[[2]], 2)
  h32_lmranks <- calculate_h_3(res[[2]], proj_model_2, 
                               n_lequal_lesser_X=n_lequal_lesser_X_2)
  expect_equivalent(h32_lmranks, H32)
})

test_that("vcov produces correct asymptotic variance estimate of rank-rank slope with no covariates", {
  load(test_path("testdata", "grouped_lmranks_cov_sigmahat_covariates_FALSE.rda"))
  res <- grouped_lmranks(r(Y) ~ r(X), data=data.frame(Y=Y,X=X), grouping_factor=G, omega=1)
  cov_mat <- vcov(res)
  sigma2hat.grouped_lmranks <- c(cov_mat[2,2]*n, cov_mat[4,4]*n)
  expect_equal(sigma2hat, sigma2hat.grouped_lmranks, tolerance=1e-5)
})

test_that("h1 works for ranked regressor with covariates and grouping", {
  load(test_path("testdata", "grouped_lmranks_cov_sigmahat_covariates_TRUE.rda"))
  res <- grouped_lmranks(r(Y) ~ r(X)+W-1, data=data.frame(Y=Y,X=X), grouping_factor=G, omega=1)
  proj_model_1 <- get_projection_model(res[[1]], 1)
  g1 <- as.integer(attr(res, "grouping_factor")) == 1
  h11_lmranks <- calculate_h_1(res[[1]], proj_model_1, g1)
  expect_equivalent(h11_lmranks, H11)
  
  proj_model_2 <- get_projection_model(res[[2]], 1)
  g2 <- as.integer(attr(res, "grouping_factor")) == 2
  h12_lmranks <- calculate_h_1(res[[2]], proj_model_2, g2)
  expect_equivalent(h12_lmranks, H12)
})

test_that("h2 works for ranked regressor with covariates and grouping", {
  load(test_path("testdata", "grouped_lmranks_cov_sigmahat_covariates_TRUE.rda"))
  res <- grouped_lmranks(r(Y) ~ r(X)+W-1, data=data.frame(Y=Y,X=X), grouping_factor=G, omega=1)
  
  n_lequal_lesser_X_1 <- count_lequal_lesser(X, X[G==1], 
                                             return_inverse_ranking=TRUE)
  n_lequal_lesser_Y_1 <- count_lequal_lesser(Y, Y[G==1], 
                                             return_inverse_ranking=TRUE)
  proj_model_1 <- get_projection_model(res[[1]], 1)
  h21_lmranks <- calculate_h_2(res[[1]], proj_model_1, 
                               n_lequal_lesser_X=n_lequal_lesser_X_1, 
                               n_lequal_lesser_Y=n_lequal_lesser_Y_1)
  expect_equivalent(h21_lmranks, H21)
  
  n_lequal_lesser_X_2 <- count_lequal_lesser(X, X[G==2], 
                                             return_inverse_ranking=TRUE)
  n_lequal_lesser_Y_2 <- count_lequal_lesser(Y, Y[G==2], 
                                             return_inverse_ranking=TRUE)
  proj_model_2 <- get_projection_model(res[[2]], 1)
  h22_lmranks <- calculate_h_2(res[[2]], proj_model_2, 
                               n_lequal_lesser_X=n_lequal_lesser_X_2, 
                               n_lequal_lesser_Y=n_lequal_lesser_Y_2)
  expect_equivalent(h22_lmranks, H22)
})

test_that("h3 works for ranked regressor with covariates and grouping", {
  load(test_path("testdata", "grouped_lmranks_cov_sigmahat_covariates_TRUE.rda"))
  res <- grouped_lmranks(r(Y) ~ r(X)+W-1, data=data.frame(Y=Y,X=X), grouping_factor=G, omega=1)
  
  n_lequal_lesser_X_1 <- count_lequal_lesser(X, X[G==1], 
                                             return_inverse_ranking=TRUE)
  proj_model_1 <- get_projection_model(res[[1]], 1)
  h31_lmranks <- calculate_h_3(res[[1]], proj_model_1, 
                               n_lequal_lesser_X=n_lequal_lesser_X_1)
  expect_equivalent(h31_lmranks, H31)
  
  n_lequal_lesser_X_2 <- count_lequal_lesser(X, X[G==2], 
                                             return_inverse_ranking=TRUE)
  proj_model_2 <- get_projection_model(res[[2]], 1)
  h32_lmranks <- calculate_h_3(res[[2]], proj_model_2, 
                               n_lequal_lesser_X=n_lequal_lesser_X_2)
  expect_equivalent(h32_lmranks, H32)
})


test_that("vcov produces correct asymptotic variance estimate of rank-rank slope with covariates", {
  load(test_path("testdata", "grouped_lmranks_cov_sigmahat_covariates_TRUE.rda"))
  res <- grouped_lmranks(r(Y) ~ r(X) + W - 1, data=data.frame(Y=Y,X=X), grouping_factor=G, omega=1)
  cov_mat <- vcov(res)
  sigma2hat.grouped_lmranks <- c(cov_mat[1,1]*n, cov_mat[5,5]*n)
  expect_equal(sigma2hat, sigma2hat.grouped_lmranks, tolerance=1e-5)
  
  W_no_intercept <- W[,-1]
  res <- grouped_lmranks(r(Y) ~ r(X) + W_no_intercept, data=data.frame(Y=Y,X=X), grouping_factor=G, omega=1)
  cov_mat <- vcov(res)
  sigma2hat.grouped_lmranks <- c(cov_mat[2,2]*n, cov_mat[6,6]*n)
  expect_equal(sigma2hat, sigma2hat.grouped_lmranks, tolerance=1e-5)
})
