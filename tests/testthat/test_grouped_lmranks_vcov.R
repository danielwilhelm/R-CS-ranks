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
  cov_mats <- calculate_grouped_lmranks_covariances(res)
  sigma2hat.grouped_lmranks <- c(cov_mats[[1]][2,2]*n, cov_mats[[2]][2,2]*n)
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


test_that("vcov produces correct asymptotic variance estimate of rank-rank slope with no covariates", {
  load(test_path("testdata", "grouped_lmranks_cov_sigmahat_covariates_TRUE.rda"))
  res <- grouped_lmranks(r(Y) ~ r(X) + W - 1, data=data.frame(Y=Y,X=X), grouping_factor=G, omega=1)
  cov_mats <- calculate_grouped_lmranks_covariances(res)
  sigma2hat.grouped_lmranks <- c(cov_mats[[1]][1,1]*n, cov_mats[[2]][1,1]*n)
  expect_equal(sigma2hat, sigma2hat.grouped_lmranks, tolerance=1e-5)
  
  W_no_intercept <- W[,-1]
  res <- grouped_lmranks(r(Y) ~ r(X) + W_no_intercept, data=data.frame(Y=Y,X=X), grouping_factor=G, omega=1)
  cov_mats <- calculate_grouped_lmranks_covariances(res)
  sigma2hat.grouped_lmranks <- c(cov_mats[[1]][2,2]*n, cov_mats[[2]][2,2]*n)
  expect_equal(sigma2hat, sigma2hat.grouped_lmranks, tolerance=1e-5)
})
