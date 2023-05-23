compare_for_tests <- function(i,j,v,omega=1){
  omega * (v[i] <= v[j]) + (1-omega) * (v[i] < v[j])
}


test_that("summary does not raise errors", {
  mod <- lmranks(r(mpg) ~ r(cyl) + disp, data=mtcars)
  expect_warning(summary(mod), "degrees of freedom")
})

test_that("vcov passes shallow checks", {
  model <- lmranks(r(mpg) ~ r(cyl) + disp, data=mtcars)
  V <- vcov(model)
  
  expect_true(isSymmetric(V))
  vals <- eigen(V, only.values = TRUE)
  expect_true(all(vals$values> 0))
  expect_equal(colnames(V),
               c("(Intercept)", "r(cyl)", "disp"))
  expect_equal(rownames(V),
               c("(Intercept)", "r(cyl)", "disp"))
})

test_that("vcov passes shallow checks in no ranked regressors case", {
  model <- lmranks(r(mpg) ~ cyl + disp, data=mtcars)
  V <- vcov(model)
  
  expect_true(isSymmetric(V))
  vals <- eigen(V, only.values = TRUE)
  expect_true(all(vals$values> 0))
})

test_that("vcov works for singular model matrix", {
  data(mtcars)
  W <- cbind(mtcars$disp, mtcars$disp)
  mod <- lmranks(r(mpg) ~ r(cyl) + W, data=mtcars)
  cov_mat <- vcov(mod)
  
  expect_true(all(!is.na(cov_mat[1:3, 1:3])))
  expect_true(all(is.na(cov_mat[4,])))
  expect_true(all(is.na(cov_mat[,4])))
})

test_that("vcov works for singular model matrix, complete=FALSE", {
  data(mtcars)
  W <- cbind(mtcars$disp, mtcars$disp)
  mod <- lmranks(r(mpg) ~ r(cyl) + W, data=mtcars)
  cov_mat <- vcov(mod, complete=FALSE)
  
  expect_true(all(!is.na(cov_mat)))
  expect_equal(c(nrow(cov_mat), ncol(cov_mat)),
               c(3,3))
})

test_that("get_and_separate_regressors works",{
  data(mtcars)
  model_1 <- lm(mpg ~ disp + cyl + hp, data=mtcars)
  model_1$rank_terms_indices <- 1
  expected_W <- cbind(rep(1, nrow(mtcars)),
                      as.matrix(mtcars[,c("cyl", "hp")]))
  colnames(expected_W) <- c("(Intercept)", "cyl", "hp")
  expected_RX <- mtcars$disp
  names(expected_RX) <- rownames(mtcars)
  expected_out_1 <- list(RX = expected_RX,
                         W = expected_W,
                         rank_column_index = 2) # Intercept
  expect_equal(get_and_separate_regressors(model_1),
               expected_out_1)
  
  W <- as.matrix(mtcars[,c("cyl", "hp")])
  RX <- mtcars$disp
  Y <- mtcars$mpg
  model_2 <- lm(Y ~ W + RX)
  model_2$rank_terms_indices <- 2
  expected_out_2 <- expected_out_1
  names(expected_out_2$RX) <- 1:length(RX)
  rownames(expected_out_2$W) <- 1:length(RX)
  colnames(expected_out_2$W) <- c("(Intercept)", "Wcyl", "Whp")
  expected_out_2$rank_column_index <- 4
  expect_equal(get_and_separate_regressors(model_2),
               expected_out_2)
})
test_that("get_and_separate_regressors works for no ranked regressors",{
  data(mtcars)
  model <- lm(mpg ~ disp + cyl + hp, data=mtcars)
  model$rank_terms_indices <- numeric(0)
  
  expected_W <- cbind(rep(1, nrow(mtcars)),
                      as.matrix(mtcars[,c("disp", "cyl", "hp")]))
  colnames(expected_W) <- c("(Intercept)", "disp", "cyl", "hp")
  expected_out <- list(RX = integer(0),
                       W = expected_W,
                       rank_column_index = integer(0))
  
  expect_equal(get_and_separate_regressors(model),
               expected_out)
})

test_that("get_projection_model works for ranked target", {
  data(mtcars)
  original_model <- lmranks(r(mpg) ~ r(disp) + cyl + hp - 1, data=mtcars,
                            omega = 1)
  
  RX <- original_model$model$`r(disp)`
  names(RX) <- rownames(mtcars)
  W <- as.matrix(mtcars[,c("cyl", "hp")])
  expected_model <- lm(RX ~ W - 1)
  expected_model$rank_terms_indices <- integer(0)
  expected_model$ranked_response <- TRUE
  expected_model$omega <- 1
  
  expect_equal(get_projection_model(original_model, 1),
               expected_model)
})

test_that("get_projection_model works for usual target with ranked regressor present", {
  data(mtcars)
  original_model <- lmranks(r(mpg) ~ r(disp) + cyl + hp - 1, data=mtcars,
                            omega = 1)
  
  RX <- original_model$model$`r(disp)`
  names(RX) <- rownames(mtcars)
  W_l <- mtcars$cyl
  names(W_l) <- rownames(mtcars)
  W_minus_l <- as.matrix(mtcars[,"hp"])
  colnames(W_minus_l) <- "hp"
  rownames(W_minus_l) <- rownames(mtcars)
  expected_model <- lm(W_l ~ RX + W_minus_l - 1)
  expected_model$rank_terms_indices <- 1
  expected_model$ranked_response <- FALSE
  expected_model$omega <- 1
  
  expect_equal(get_projection_model(original_model, 2),
               expected_model)
})

test_that("get_projection_model works for usual target with ranked regressor present", {
  data(mtcars)
  original_model <- lmranks(r(mpg) ~ r(disp) + cyl + hp - 1, data=mtcars,
                            omega = 1)
  
  RX <- original_model$model$`r(disp)`
  names(RX) <- rownames(mtcars)
  W_l <- mtcars$cyl
  names(W_l) <- rownames(mtcars)
  W_minus_l <- as.matrix(mtcars[,"hp"])
  colnames(W_minus_l) <- "hp"
  rownames(W_minus_l) <- rownames(mtcars)
  expected_model <- lm(W_l ~ RX + W_minus_l - 1)
  expected_model$rank_terms_indices <- 1
  expected_model$ranked_response <- FALSE
  expected_model$omega <- 1
  
  expect_equal(get_projection_model(original_model, 2),
               expected_model)
})

test_that("get_projection_model works with no ranked regressors", {
  data(mtcars)
  original_model <- lmranks(r(mpg) ~ disp + cyl + hp - 1, data=mtcars,
                            omega=1)
  
  W_l <- mtcars$cyl
  names(W_l) <- rownames(mtcars)
  W_minus_l <- as.matrix(mtcars[,c("disp", "hp")])
  colnames(W_minus_l) <- c("disp", "hp")
  rownames(W_minus_l) <- rownames(mtcars)
  expected_model <- lm(W_l ~ W_minus_l - 1)
  expected_model$rank_terms_indices <- integer(0)
  expected_model$ranked_response <- FALSE
  expected_model$omega <- 1
  
  expect_equal(get_projection_model(original_model, 2),
               expected_model)
})

test_that("get_projection_model works for intercept", {
  data(mtcars)
  original_model <- lmranks(r(mpg) ~ r(disp) + cyl + hp, data=mtcars,
                            omega=1)
  
  RX <- original_model$model$`r(disp)`
  names(RX) <- rownames(mtcars)
  W_l <- rep(1, nrow(mtcars))
  names(W_l) <- rownames(mtcars)
  W_minus_l <- as.matrix(mtcars[,c("cyl", "hp")])
  expected_model <- lm(W_l ~ RX + W_minus_l - 1)
  expected_model$rank_terms_indices <- 1
  expected_model$ranked_response <- FALSE
  expected_model$omega <- 1
  
  expect_equal(get_projection_model(original_model, 1),
               expected_model)
})

test_that("get_projection_model works for model with 1 usual regressor", {
  data(mtcars)
  original_model <- lmranks(r(mpg) ~ r(disp), data=mtcars,
                            omega=1)
  
  RX <- original_model$model$`r(disp)`
  names(RX) <- rownames(mtcars)
  W_l <- rep(1, nrow(mtcars))
  names(W_l) <- rownames(mtcars)
  expected_model <- lm(W_l ~ RX - 1)
  expected_model$rank_terms_indices <- 1
  expected_model$ranked_response <- FALSE
  expected_model$omega <- 1
  
  expect_equal(get_projection_model(original_model, 1),
               expected_model)
})

test_that("calculate_g_l_3 works for no ranked regressors case", {
  data(mtcars)
  n_lequal_lesser_X <- NULL
  W <- as.matrix(mtcars[, c("cyl", "hp")])
  W <- cbind(rep(1, nrow(W)),
             W)
  
  original_model <- lmranks(r(mpg) ~ cyl + hp, data=mtcars)
  RY <- original_model$model$`r(mpg)`
  W_l <- W[,2]
  W_minus_l <- W[,-2]
  proj_model <- lm(W_l ~ W_minus_l - 1)
  proj_model$rank_terms_indices <- integer(0)
  proj_model$ranked_response <- FALSE
  
  expected_out <- rep(0, nrow(mtcars))
  
  expect_equal(calculate_g_l_3(original_model, proj_model, n_lequal_lesser_X=n_lequal_lesser_X),
               expected_out)
})

test_that("calculate_g_l_3 works for proj_model with ranked response", {
  data(mtcars)
  original_model <- lmranks(r(mpg) ~ r(cyl) + hp, data=mtcars, omega=1)
  RX <- original_model$model$`r(cyl)`
  n_lequal_lesser_X <- count_lequal_lesser(mtcars$cyl, return_inverse_ranking = TRUE)
  expected_I_X <- function(i,j) compare_for_tests(i,j,mtcars$cyl)
  W <- as.matrix(mtcars[, c("hp")])
  W <- cbind(rep(1, nrow(W)),
             W)
  
  proj_model <- lm(RX ~ W - 1)
  proj_model$rank_terms_indices <- integer(0)
  proj_model$ranked_response <- TRUE
  proj_model$omega <- 1
  
  original_resid <- resid(original_model)
  predictor <- proj_model$fitted.values
  
  expected_out <- sapply(1:length(RX), function(i){
    mean(
      sapply(1:length(RX), function(j){
        original_resid[j] * (expected_I_X(i,j) - predictor[j])
      })
    )
  })
  
  expect_equal(calculate_g_l_3(original_model, proj_model, n_lequal_lesser_X=n_lequal_lesser_X),
               expected_out)
})

test_that("calculate_g_l_3 works for proj_model with usual response", {
  data(mtcars)
  original_model <- lmranks(r(mpg) ~ r(cyl) + hp, data=mtcars, omega=1)
  RX <- original_model$model$`r(cyl)`
  n_lequal_lesser_X <- count_lequal_lesser(mtcars$cyl, return_inverse_ranking = TRUE)
  expected_I_X <- function(i,j) compare_for_tests(i,j,mtcars$cyl)
  W_l <- mtcars$hp
  W_minus_l <- rep(1, length(W_l))
  
  
  proj_model <- lm(W_l ~ RX + W_minus_l - 1)
  proj_model$rank_terms_indices <- 1
  proj_model$ranked_response <- FALSE
  proj_model$omega <- 1
  
  original_resid <- resid(original_model)
  coefs <- coef(proj_model)
  proj_usual_predictor <- W_minus_l * coefs[2]
  
  expected_out <- sapply(1:length(RX), function(i){
    mean(
      sapply(1:length(RX), function(j){
        original_resid[j] * (W_l[j] - coefs[1]*expected_I_X(i,j) - proj_usual_predictor[j])
      })
    )
  })
  
  expect_equal(calculate_g_l_3(original_model, proj_model, n_lequal_lesser_X=n_lequal_lesser_X),
               expected_out)
})

test_that("calculate_weighted_ineq_resid_means", {
  skip()
})

test_that("extract_nonrank_predictor works when no ranked regressors present", {
  data(mtcars)
  model <- lmranks(r(mpg) ~ cyl + disp, data=mtcars)
  expected <- model$fitted.values
  names(expected) <- NULL
  expect_equal(extract_nonrank_predictor(model),
               expected)
})

test_that("extract_nonrank_predictor works in singular case", {
  data(mtcars)
  mtcars$disp2 <- mtcars$disp
  model <- lmranks(r(mpg) ~ disp + disp2, data=mtcars)
  expected_nonrank_predictor <- coef(model)[1] + coef(model)[2] * mtcars$disp
  expect_equal(extract_nonrank_predictor(model),
               expected_nonrank_predictor)
})

test_that("extract_nonrank_predictor works when ranked regressors present", {
  data(mtcars)
  model <- lmranks(r(mpg) ~ r(cyl) + disp, data=mtcars)
  expected_nonrank_predictor <- coef(model)[1] + coef(model)[3] * mtcars$disp
  expect_equal(extract_nonrank_predictor(model),
               expected_nonrank_predictor)
})

test_that("get_ineq_indicator_function returns functions with identical arguments", {
  v <- 1:10
  f1 <- get_ineq_indicator_function(FALSE, v)
  f2 <- get_ineq_indicator_function(TRUE, v)
  
  expect_equal(args(f1), args(f2))
})

test_that("get_ineq_indicator_function remembers vector values across calls", {
  v1 <- 1:3
  v2 <- 4:6
  v3 <- 7:9
  
  f1 <- get_ineq_indicator_function(FALSE, v1)
  f2 <- get_ineq_indicator_function(FALSE, v2)
  f3 <- get_ineq_indicator_function(FALSE, v3)
  
  expect_equal(f1(1,2,3), v1)
  expect_equal(f2(1,2,3), v2)
  expect_equal(f2(1,2,3), v2)
})

test_that("get_ineq_indicator works for sorted data", {
  original_v <- c(1,3,4,4,4,7,7,10)
  n_lequal_lesser <- matrix(c(
    1,0,1,
    2,1,2,
    5,2,3,
    5,2,4,
    5,2,5,
    7,5,6,
    7,5,7,
    8,7,8
  ), byrow=TRUE, ncol = 3)
  for(i in 1:length(original_v)){
    expected_om0 <- sapply(1:length(original_v), function(j)
      compare_for_tests(i,j,original_v,omega=0))
    expected_om0.4 <- sapply(1:length(original_v), function(j)
      compare_for_tests(i,j,original_v,omega=0.4))
    expected_om1 <- sapply(1:length(original_v), function(j)
      compare_for_tests(i,j,original_v,omega=1))
    
    expect_equal(get_ineq_indicator(n_lequal_lesser, i, 0),
                 expected_om0)
    expect_equal(get_ineq_indicator(n_lequal_lesser, i, 0.4),
                 expected_om0.4)
    expect_equal(get_ineq_indicator(n_lequal_lesser, i, 1),
                 expected_om1)
  }
})

test_that("get_ineq_indicator works for unsorted data", {
  original_v <- c(4,4,4,3,1,10,7,7)
  n_lequal_lesser <- matrix(c(
    5,2,3,
    5,2,4,
    5,2,5,
    2,1,2,
    1,0,1,
    8,7,8,
    7,5,6,
    7,5,7
  ), byrow=TRUE, ncol = 3)
  for(i in 1:length(original_v)){
    expected_om0 <- sapply(1:length(original_v), function(j)
      compare_for_tests(i,j,original_v,omega=0))
    expected_om0.4 <- sapply(1:length(original_v), function(j)
      compare_for_tests(i,j,original_v,omega=0.4))
    expected_om1 <- sapply(1:length(original_v), function(j)
      compare_for_tests(i,j,original_v,omega=1))
    
    expect_equal(get_ineq_indicator(n_lequal_lesser, i, 0),
                 expected_om0)
    expect_equal(get_ineq_indicator(n_lequal_lesser, i, 0.4),
                 expected_om0.4)
    expect_equal(get_ineq_indicator(n_lequal_lesser, i, 1),
                 expected_om1)
  }
})

test_that("replace_ranks_with_ineq_indicator_and_calculate_residuals works 
          for usual regressors", {
            data(mtcars)
            i <- 5
            RY <- frank(mtcars$mpg)
            I_Y <- sapply(1:length(mtcars$mpg), function(j)
              compare_for_tests(i, j, mtcars$mpg))
            W <- as.matrix(mtcars[, c("disp", "hp")])
            W <- cbind(rep(1, nrow(W)),
                       W)
            original_model <- lm(RY ~ W - 1)
            original_model$rank_terms_indices <- integer(0)
            original_model$ranked_response <- TRUE
            
            coefs <- coef(original_model)
            predictor <- W %*% coefs
            expected_resid <- I_Y - predictor
            
            expect_equal(replace_ranks_with_ineq_indicator_and_calculate_residuals(
              original_model, nonrank_predictor = predictor, I_Y = I_Y
            ), expected_resid)
            
          })

test_that("replace_ranks_with_ineq_indicator_and_calculate_residuals works 
          for ranked regressors", {
            data(mtcars)
            i <- 5
            RX <- frank(mtcars$cyl)
            I_X <- sapply(1:length(mtcars$cyl), function(j)
              compare_for_tests(i, j, mtcars$cyl))
            RY <- frank(mtcars$mpg)
            I_Y <- sapply(1:length(mtcars$mpg), function(j)
              compare_for_tests(i, j, mtcars$mpg))
            W <- as.matrix(mtcars[, c("hp")])
            W <- cbind(rep(1, nrow(W)),
                       W)
            original_model <- lm(RY ~ RX + W - 1)
            original_model$rank_terms_indices <- 1
            original_model$ranked_response <- TRUE
            
            coefs <- coef(original_model)
            predictor <- W %*% coefs[-1]
            expected_resid <- I_Y - coefs[1] * I_X - predictor
            
            expect_equal(replace_ranks_with_ineq_indicator_and_calculate_residuals(
              original_model, nonrank_predictor = predictor, I_X = I_X, I_Y = I_Y
            ), expected_resid)
          })
test_that("vcov produces correct asymptotic variance estimate of rank-rank slope with covariates present", {
  load(test_path("testdata", "lmranks_cov_sigmahat_covariates_TRUE.rda"))
  res <- lmranks(r(Y) ~ r(X) + W - 1)
  sigma2hat.lmranks <- vcov(res)[1,1]*n
  expect_equal(sigma2hat, sigma2hat.lmranks, tolerance=1e-5)
})

test_that("vcov produces correct asymptotic variance estimate of rank-rank slope with no covariates", {
  load(test_path("testdata", "lmranks_cov_sigmahat_covariates_FALSE.rda"))
  res <- lmranks(r(Y) ~ r(X))
  sigma2hat.lmranks <- vcov(res)[2,2]*n
  expect_equal(sigma2hat, sigma2hat.lmranks, tolerance=1e-5)
})

test_that("vcov produces asymptotic variance estimate of rank-rank slope equal to that of Hoeffding (1948)", {
  load(test_path("testdata", "lmranks_cov_sigmahat_Hoefding.rda"))
  res <- lmranks(r(Y) ~ r(X))
  sigma2hat.lmranks <- vcov(res)[2,2]*n
  expect_equal(sigma2hat, sigma2hat.lmranks, tolerance=1e-3)
})



