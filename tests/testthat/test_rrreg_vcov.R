test_that("summary does not raise errors", {
  mod <- lmranks(r(mpg) ~ r(cyl) + disp, data=mtcars)
  expect_silent(summary(mod))
})

test_that("vcov passes shallow checks", {
  model <- lmranks(r(mpg) ~ r(cyl) + disp, data=mtcars)
  V <- vcov(model)
  
  expect_true(isSymmetric(V))
  vals <- eigen(V, only.values = TRUE)
  expect_true(all(vals > 0))
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
  original_model <- lmranks(r(mpg) ~ r(disp) + cyl + hp - 1, data=mtcars)
  
  RX <- frank(mtcars$disp)
  names(RX) <- rownames(mtcars)
  W <- as.matrix(mtcars[,c("cyl", "hp")])
  expected_model <- lm(RX ~ W - 1)
  expected_model$rank_terms_indices <- integer(0)
  expected_model$ranked_response <- TRUE
  
  expect_equal(get_projection_model(original_model, 1),
               expected_model)
})

test_that("get_projection_model works for usual target", {
  data(mtcars)
  original_model <- lmranks(r(mpg) ~ r(disp) + cyl + hp - 1, data=mtcars)
  
  RX <- frank(mtcars$disp)
  names(RX) <- rownames(mtcars)
  W_l <- mtcars$cyl
  names(W_l) <- rownames(mtcars)
  W_minus_l <- as.matrix(mtcars[,"hp"])
  colnames(W_minus_l) <- "hp"
  rownames(W_minus_l) <- rownames(mtcars)
  expected_model <- lm(W_l ~ RX + W_minus_l - 1)
  expected_model$rank_terms_indices <- 1
  expected_model$ranked_response <- FALSE
  
  expect_equal(get_projection_model(original_model, 2),
               expected_model)
})

test_that("get_projection_model works for intercept", {
  data(mtcars)
  original_model <- lmranks(r(mpg) ~ r(disp) + cyl + hp, data=mtcars)
  
  RX <- frank(mtcars$disp)
  names(RX) <- rownames(mtcars)
  W_l <- rep(1, nrow(mtcars))
  names(W_l) <- rownames(mtcars)
  W_minus_l <- as.matrix(mtcars[,c("cyl", "hp")])
  expected_model <- lm(W_l ~ RX + W_minus_l - 1)
  expected_model$rank_terms_indices <- 1
  expected_model$ranked_response <- FALSE
  
  expect_equal(get_projection_model(original_model, 1),
               expected_model)
})

test_that("get_projection_model works for model with 1 usual regressor", {
  data(mtcars)
  original_model <- lmranks(r(mpg) ~ r(disp), data=mtcars)
  
  RX <- frank(mtcars$disp)
  names(RX) <- rownames(mtcars)
  W_l <- rep(1, nrow(mtcars))
  names(W_l) <- rownames(mtcars)
  expected_model <- lm(W_l ~ RX - 1)
  expected_model$rank_terms_indices <- 1
  expected_model$ranked_response <- FALSE
  
  expect_equal(get_projection_model(original_model, 1),
               expected_model)
})

test_that("calculate_g_l_3 works for proj_model with ranked response", {
  data(mtcars)
  RX <- frank(mtcars$cyl)
  I_X <- compare(mtcars$cyl)
  RY <- frank(mtcars$mpg)
  W <- as.matrix(mtcars[, c("hp")])
  W <- cbind(rep(1, nrow(W)),
             W)
  
  original_model <- lmranks(r(mpg) ~ r(cyl) + hp, data=mtcars)
  proj_model <- lm(RX ~ W - 1)
  proj_model$rank_terms_indices <- integer(0)
  proj_model$ranked_response <- TRUE
  
  original_resid <- resid(original_model)
  predictor <- proj_model$fitted.values
  
  expected_out <- sapply(1:length(RY), function(i){
    mean(
      sapply(1:length(RY), function(j){
        original_resid[j] * (I_X[i,j] - predictor[j])
      })
    )
  })
  
  expect_equal(calculate_g_l_3(original_model, proj_model, I_X=I_X),
               expected_out)
})

test_that("calculate_g_l_3 works for proj_model with usual response", {
  data(mtcars)
  RX <- frank(mtcars$cyl)
  I_X <- compare(mtcars$cyl)
  RY <- frank(mtcars$mpg)
  I_Y <- compare(mtcars$mpg)
  W_l <- mtcars$hp
  W_minus_l <- rep(1, length(W_l))
  
  original_model <- lmranks(r(mpg) ~ r(cyl) + hp, data=mtcars)
  proj_model <- lm(W_l ~ RX + W_minus_l - 1)
  proj_model$rank_terms_indices <- 1
  proj_model$ranked_response <- FALSE
  
  original_resid <- resid(original_model)
  coefs <- coef(proj_model)
  proj_usual_predictor <- W_minus_l * coefs[2]
  
  expected_out <- sapply(1:length(RY), function(i){
    mean(
      sapply(1:length(RY), function(j){
        original_resid[j] * (W_l[j] - coefs[1]*I_X[i,j] - proj_usual_predictor[j])
      })
    )
  })
  
  expect_equal(calculate_g_l_3(original_model, proj_model, I_X=I_X),
               expected_out)
})

test_that("replace_ranks_with_ineq_indicator_and_calculate_residuals works 
          for ranked response, usual regressors", {
  data(mtcars)
  RY <- frank(mtcars$mpg)
  I_Y <- compare(mtcars$mpg)
  W <- as.matrix(mtcars[, c("disp", "hp")])
  W <- cbind(rep(1, nrow(W)),
             W)
  original_model <- lm(RY ~ W - 1)
  original_model$rank_terms_indices <- integer(0)
  original_model$ranked_response <- TRUE
  
  coefs <- coef(original_model)
  predictor <- W %*% coefs
  expected_resid <- sapply(1:length(RY), function(j){
    I_Y[,j] - predictor[j]
  })
  expected_resid <- t(expected_resid)
  
  expect_equal(replace_ranks_with_ineq_indicator_and_calculate_residuals(
    original_model, I_Y = I_Y
  ), expected_resid)
  
})

test_that("replace_ranks_with_ineq_indicator_and_calculate_residuals works 
          for ranked response, usual regressors", {
  data(mtcars)
  RY <- frank(mtcars$mpg)
  I_Y <- compare(mtcars$mpg)
  W <- as.matrix(mtcars[, c("disp", "hp")])
  W <- cbind(rep(1, nrow(W)),
             W)
  original_model <- lm(RY ~ W - 1)
  original_model$rank_terms_indices <- integer(0)
  original_model$ranked_response <- TRUE
  
  coefs <- coef(original_model)
  predictor <- W %*% coefs
  expected_resid <- sapply(1:length(RY), function(j){
    I_Y[,j] - predictor[j]
  })
  expected_resid <- t(expected_resid)
  
  expect_equal(replace_ranks_with_ineq_indicator_and_calculate_residuals(
    original_model, I_Y = I_Y
  ), expected_resid)
  
})

test_that("replace_ranks_with_ineq_indicator_and_calculate_residuals works 
          for usual response, ranked regressors", {
    data(mtcars)
    RX <- frank(mtcars$cyl)
    I_X <- compare(mtcars$cyl)
    Y <- frank(mtcars$mpg)
    W <- as.matrix(mtcars[, c("hp")])
    W <- cbind(rep(1, nrow(W)),
               W)
    original_model <- lm(Y ~ RX + W - 1)
    original_model$rank_terms_indices <- 1
    original_model$ranked_response <- FALSE
    
    coefs <- coef(original_model)
    predictor <- W %*% coefs[-1]
    expected_resid <- sapply(1:length(Y), function(j){
      Y[j] - coefs[1] * I_X[,j] - predictor[j]
    })
    expected_resid <- t(expected_resid)
    
    expect_equal(replace_ranks_with_ineq_indicator_and_calculate_residuals(
      original_model, I_X = I_X
    ), expected_resid)
})

test_that("replace_ranks_with_ineq_indicator_and_calculate_residuals works 
          for usual response, ranked regressors", {
    data(mtcars)
    RX <- frank(mtcars$cyl)
    I_X <- compare(mtcars$cyl)
    RY <- frank(mtcars$mpg)
    I_Y <- compare(mtcars$mpg)
    W <- as.matrix(mtcars[, c("hp")])
    W <- cbind(rep(1, nrow(W)),
               W)
    original_model <- lm(RY ~ RX + W - 1)
    original_model$rank_terms_indices <- 1
    original_model$ranked_response <- TRUE
    
    coefs <- coef(original_model)
    predictor <- W %*% coefs[-1]
    expected_resid <- sapply(1:length(RY), function(j){
      I_Y[,j] - coefs[1] * I_X[,j] - predictor[j]
    })
    expected_resid <- t(expected_resid)
    
    expect_equal(replace_ranks_with_ineq_indicator_and_calculate_residuals(
      original_model, I_X = I_X, I_Y = I_Y
    ), expected_resid)
})