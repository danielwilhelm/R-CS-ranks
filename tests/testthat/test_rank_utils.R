### compare ###

test_that("compare works for sorted input", {
  x <- c(1,3,4,4,4,7,7,10)
  expected_compare_om0.4 <- matrix(c(
    0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
    1.0, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
    1.0, 1.0, 0.4, 0.4, 0.4, 0.0, 0.0, 0.0,
    1.0, 1.0, 0.4, 0.4, 0.4, 0.0, 0.0, 0.0,
    1.0, 1.0, 0.4, 0.4, 0.4, 0.0, 0.0, 0.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 0.4, 0.4, 0.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 0.4, 0.4, 0.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.4
  ), byrow = TRUE, ncol = 8)
  expected_compare_om0 <- expected_compare_om0.4
  expected_compare_om0[expected_compare_om0 == 0.4] <- 0
  expected_compare_om1 <- expected_compare_om0.4
  expected_compare_om1[expected_compare_om1 == 0.4] <- 1
  
  expect_equal(compare(x, omega = 0, increasing = TRUE, na.rm = FALSE),
               expected_compare_om0)
  expect_equal(compare(x, omega = 0.4, increasing = TRUE, na.rm = FALSE),
               expected_compare_om0.4)
  expect_equal(compare(x, omega = 1, increasing = TRUE, na.rm = FALSE),
               expected_compare_om1)
})

test_that("compare works for unsorted input", {
  x <- c(4, 3, 4, 4, 7, 7, 10, 1)
  expected_compare_om0.4 <- matrix(c(
    0.4, 1.0, 0.4, 0.4, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.4, 1.0, 0.4, 0.4, 0.0, 0.0, 0.0, 1.0,
    0.4, 1.0, 0.4, 0.4, 0.0, 0.0, 0.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 0.4, 0.4, 0.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 0.4, 0.4, 0.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.4, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4
  ), byrow = TRUE, ncol = 8)
  expected_compare_om0 <- expected_compare_om0.4
  expected_compare_om0[expected_compare_om0 == 0.4] <- 0
  expected_compare_om1 <- expected_compare_om0.4
  expected_compare_om1[expected_compare_om1 == 0.4] <- 1
  
  expect_equal(compare(x, omega = 0, increasing = TRUE, na.rm = FALSE),
               expected_compare_om0)
  expect_equal(compare(x, omega = 0.4, increasing = TRUE, na.rm = FALSE),
               expected_compare_om0.4)
  expect_equal(compare(x, omega = 1, increasing = TRUE, na.rm = FALSE),
               expected_compare_om1)
})

test_that("compare handles NAs", {
  x <- c(4, 3, 4, NA, 7, 7, NA, 1)
  expected_compare_na.rmT <- matrix(c(
    0.4, 1.0, 0.4, 0.0, 0.0, 1.0,
    0.0, 0.4, 0.0, 0.0, 0.0, 1.0,
    0.4, 1.0, 0.4, 0.0, 0.0, 1.0,
    1.0, 1.0, 1.0, 0.4, 0.4, 1.0,
    1.0, 1.0, 1.0, 0.4, 0.4, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.4
  ), byrow = TRUE, ncol = 6)
  expected_compare_na.rmF <- matrix(c(
    0.4, 1.0, 0.4, NA, 0.0, 0.0, NA, 1.0,
    0.0, 0.4, 0.0, NA, 0.0, 0.0, NA, 1.0,
    0.4, 1.0, 0.4, NA, 0.0, 0.0, NA, 1.0,
     NA,  NA,  NA, NA,  NA,  NA, NA,  NA, 
    1.0, 1.0, 1.0, NA, 0.4, 0.4, NA, 1.0,
    1.0, 1.0, 1.0, NA, 0.4, 0.4, NA, 1.0,
     NA,  NA,  NA, NA,  NA,  NA, NA,  NA,
    0.0, 0.0, 0.0, NA, 0.0, 0.0, NA, 0.4
  ), byrow = TRUE, ncol = 8)
  expect_equal(compare(x, omega = 0.4, increasing = TRUE, na.rm = TRUE),
               expected_compare_na.rmT)
  expect_error(compare(x, omega = 0.4, increasing = TRUE, na.rm = FALSE))
})

### irank ###

test_that("irank works for sorted input", {
  x <- c(1,3,4,4,4,7,7,10)
  expected_irank_om0 <- c(8,7,4,4,4,2,2,1)
  expected_irank_om0_increasing <- c(1,2,3,3,3,6,6,8)
  expected_irank_om1 <- c(8,7,6,6,6,3,3,1)
  expected_irank_om0.5 <- c(8,7,5,5,5,2.5,2.5,1)
  
  expect_equal(irank(x, omega = 0),
               expected_irank_om0)
  expect_equal(irank(x, omega = 1),
               expected_irank_om1)
  expect_equal(irank(x, omega = 0.5),
               expected_irank_om0.5)
  
  expect_equal(irank(x, increasing = TRUE),
               expected_irank_om0_increasing)
})

test_that("irank works for unsorted input", {
  x <- c(1,3,4,4,4,7,7,10)
  shuffle <- sample(length(x))
  x <- x[shuffle]
  expected_irank_om0 <- c(8,7,4,4,4,2,2,1)[shuffle]
  expected_irank_om0.5 <- c(8,7,5,5,5,2.5,2.5,1)[shuffle]
  expect_equal(irank(x, omega = 0.5),
               expected_irank_om0.5)
})
