test_that("vcov produces correct asymptotic variance estimate of rank-rank slope with no covariates", {
  load(test_path("testdata", "grouped_lmranks_cov_sigmahat_covariates_FALSE.rda"))
  res <- grouped_lmranks(r(Y) ~ r(X), data=data.frame(Y=Y,X=X), grouping_factor=G, omega=1)
  cov_mats <- calculate_grouped_lmranks_covariances(res)
  sigma2hat.grouped_lmranks <- c(cov_mats[[1]][2,2]*n, cov_mats[[2]][2,2]*n)
  expect_equal(sigma2hat, sigma2hat.grouped_lmranks, tolerance=1e-5)
})

test_that("vcov produces correct asymptotic variance estimate of rank-rank slope with no covariates", {
  load(test_path("testdata", "grouped_lmranks_cov_sigmahat_covariates_TRUE.rda"))
  res <- grouped_lmranks(r(Y) ~ r(X) + W - 1, data=data.frame(Y=Y,X=X), grouping_factor=G, omega=1)
  cov_mats <- calculate_grouped_lmranks_covariances(res)
  sigma2hat.grouped_lmranks <- c(cov_mats[[1]][1,1]*n, cov_mats[[2]][1,1]*n)
  expect_equal(sigma2hat, sigma2hat.grouped_lmranks, tolerance=1e-5)
})
