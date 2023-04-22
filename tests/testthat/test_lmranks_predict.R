test_that("predict works for old data supplied implicitly", {
  Y <- c(3,1,2,4,5)
  y_frank <- c(0.6, 1.0, 0.8, 0.4, 0.2)
  X <- 1:5
  omega <- 0.5
  x_frank <- c(1.0, 0.8, 0.6, 0.4, 0.2)
  W <- c(1,3,2,5,4)
  
  data <- data.frame(y = Y,
                     x = X,
                     w = W)
  model <- lmranks(r(y) ~ r(x) + w, data=data)
  expect_equal(predict(model),
               model$fitted.values)
})

test_that("predict works for old data supplied explicitly", {
  Y <- c(3,1,2,4,5)
  y_frank <- c(0.6, 1.0, 0.8, 0.4, 0.2)
  X <- 1:5
  omega <- 0.5
  x_frank <- c(1.0, 0.8, 0.6, 0.4, 0.2)
  W <- c(1,3,2,5,4)
  
  data <- data.frame(y = Y,
                     x = X,
                     w = W)
  model <- lmranks(r(y) ~ r(x) + w, data=data)
  expect_equal(predict(model, data),
               model$fitted.values)
})

test_that("predict works for new data", {
  Y <- c(3, 1, 2, 4, 5)
  X <- 1:5
  omega <- 0.5
  W <- c(1, 3, 2, 5, 4)
  
  data <- data.frame(y = Y,
                     x = X,
                     w = W)
  model <- lmranks(r(y) ~ r(x) + w, data = data, omega=omega)
  coefs <- coef(model)
  
  new_data <- data.frame(x = c(1.5, 2.5, 3),
                         w = c(1.5, 2.5, 3))
  new_x_rank <- c(0.2, 0.4, 0.5)
  expected_prediction <-
    coefs[1] + coefs[2] * new_x_rank + coefs[3] * new_data$w
  names(expected_prediction) <- 1:3
  expect_equal(predict(model, new_data),
               expected_prediction)
})

test_that("predict works for complicated response", {
  model <- lmranks(r(log(mpg)) ~ r(I(disp ^ 2)) + cyl + hp, data = mtcars)
  expect_equal(predict(model, mtcars),
               model$fitted.values)
})