# we want a very simple example

test_that("simple_lmranks works for full correlation", {
  Y <- c(1,2,3)
  X <- c(4,5,6)
  expected_rho <- 1
  expected_se <- 0
  actual <- simple_lmranks(Y, X)
  expect_equal(actual$rhohat, expected_rho)
  expect_equal(actual$se, expected_se)
})

test_that("simple_lmranks works for no correlation",{
  Y <- c(2,1,0,1,2)
  X <- 1:5
  expected_rho <- 0
  actual <- simple_lmranks(Y, X)
  expect_equal(actual$rhohat, expected_rho)
})

test_that("simple_lmranks works for NULL W", {
  Y <- c(2,0,0,3,3)
  omega <- 0.5
  y_frank <- c(0.6, 0.9, 0.9, 0.3,0.3)
  X <- 1:5
  x_frank <- (5:1) / 5
  m <- lm(y_frank ~ x_frank)
  expected_rho <- coef(m)[2]
  names(expected_rho) <- NULL
  sigma_nu <- var(x_frank)
  gammahat <- 0.6 # mean of frank of x
  betahat <- coef(m)[1]
  h1 <- (y_frank - expected_rho * x_frank - betahat)*(x_frank - gammahat)
  h2 <- sapply(1:5, function(i)
    mean(sapply(1:5, function(j){
      (omega * (Y[i] <= Y[j]) + (1-omega)*(Y[i] < Y[j]) - expected_rho * 
         (omega * (X[i] <= X[j]) + (1-omega)*(X[i] < X[j])) - betahat) * 
        (x_frank[j] - gammahat)
    })))
  h3 <- sapply(1:5, function(i)
    mean(sapply(1:5, function(j){
      (y_frank[j] - expected_rho * x_frank[j] - betahat) * 
         (omega * (X[i] <= X[j]) + (1-omega)*(X[i] < X[j]) - gammahat)
    })))
  expected_se <- sqrt(mean((h1 + h2 + h3)^2) / sigma_nu^2)
  
  actual <- simple_lmranks(Y, X, omega = omega)
  expect_equal(actual$rhohat, expected_rho)
  expect_equal(actual$se * sqrt(5), expected_se)
})