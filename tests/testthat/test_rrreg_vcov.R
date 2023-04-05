test_that("summary does not raise errors", {
  mod <- lmranks(r(mpg) ~ r(cyl) + disp, data=mtcars)
  expect_silent(summary(mod))
})

test_that("vcov works", {
  
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

test_that("calculate_g_l_3 works", {
  
})

test_that("replace_ranks_with_ineq_indicator_and_calculate_residuals works", {
  
})