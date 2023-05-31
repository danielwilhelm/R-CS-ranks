### count_lequal_lesser ###
test_that("count_lequal_lesser works for sorted x", {
  v <- c(1,3,4,4,4,7,7,10)
  x <- c(0,1,2,3,4,5,7,8,10,11)
  expected <- list(n_lequal = c(0,1,1,2,5,5,7,7,8,8),
                   n_lesser = c(0,0,1,1,2,5,5,7,7,8))
  expect_equal(count_lequal_lesser(x, v),
               expected)
  v <- c(4,4,4,3,1,10,7,7)
  expect_equal(count_lequal_lesser(x, v),
               expected)
})

test_that("count_lequal_lesser works for unsorted x", {
  v <- c(1,3,4,4,4,7,7,10)
  x <- c(5,10,8,1,7,2,3,4,11,0)
  expected <- list(n_lequal = c(5,8,7,1,7,1,2,5,8,0),
                   n_lesser = c(5,7,7,0,5,1,1,2,8,0))
  expect_equal(count_lequal_lesser(x, v),
               expected)
  v <- c(4,4,4,3,1,10,7,7)
  expect_equal(count_lequal_lesser(x, v),
               expected)
})

test_that("count_lequal_lesser with default v works", {
  x <- c(1,3,4,4,4,7,7,10)
  expect_equal(count_lequal_lesser(x, x),
               count_lequal_lesser(x))
})

test_that("count_lequal_lesser return_inverse_ranking argument works", {
  x <- c(4,4,4,3,1,10,7,7)
  x_sorted <- sort(x)
  out <- count_lequal_lesser(x, return_inverse_ranking = TRUE)
  expect_false(is.null(out[["inverse_ranking"]]))
  expect_equal(x_sorted[out$inverse_ranking], x)
})

test_that("count_lequal_lesser return_inverse_ranking argument works with custom v", {
  v <- c(4,4,4,3,1,10,7,7)
  v_sorted <- sort(v)
  x <- c(1,3,5,7,9)
  out <- count_lequal_lesser(x, v, return_inverse_ranking = TRUE)
  expect_false(is.null(out[["inverse_ranking"]]))
  expect_equal(v_sorted[out$inverse_ranking], v)
})

### compare ###
test_that("frank_against returns error for matrix input", {
  # This behavior ensures, that user cannot pass r(MATRIX) in lmranks
  expect_error(frank(matrix(1:12, ncol = 3),
                                    omega=0.4, increasing=TRUE,
                                    na.rm = FALSE))
  expect_error(process_compare_args(matrix(1:12, ncol = 3),
                                    omega=0.4, increasing=TRUE,
                                    na.rm = FALSE))
})

### irank ###

test_that("irank works for sorted input", {
  x <- c(1,3,4,4,4,7,7,10)
  expected_irank_om0 <- c(8,7,4,4,4,2,2,1)
  expected_irank_om1 <- c(8,7,6,6,6,3,3,1)
  expected_irank_om0.5 <- c(8,7,5,5,5,2.5,2.5,1)
  
  expect_equal(irank(x, omega = 0),
               expected_irank_om0)
  expect_equal(irank(x, omega = 1),
               expected_irank_om1)
  expect_equal(irank(x, omega = 0.5),
               expected_irank_om0.5)
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

test_that("irank's increasing argument works", {
  x <- c(1,3,4,4,4,7,7,10)
  x_neg <- -x
  expected_irank_om0.4 <- c(1,2,3.8,3.8,3.8,6.4,6.4,8)
  
  expect_equal(irank(x, increasing = TRUE, omega=0.4),
               expected_irank_om0.4)
  expect_equal(irank(x_neg, increasing = FALSE, omega=0.4),
               expected_irank_om0.4)
})

test_that("irank_against works", {
  v <- c(1,3,4,4,4,7,7,10)
  x <- c(0,1,2,3,4,5,7,8,10,11)
  expected_irank <- c(8.5, 8, 7.5, 7, 5, 3.5, 2.5, 1.5, 1, 0.5)
  expect_equal(irank_against(x, v, omega = 0.5),
               expected_irank)
})

test_that("irank_against's increasing argument works", {
  v <- c(1,3,4,4,4,7,7,10)
  x <- c(0,1,2,3,4,5,7,8,10,11)
  expected_irank <- c(0.5, 1, 1.5, 2, 4, 5.5, 6.5, 7.5, 8, 8.5)
  expect_equal(irank_against(x, v, omega = 0.5, increasing = TRUE),
               expected_irank)
})

test_that("irank_against handles NAs", {
  v <- c(1,3,4,4,4,7,7,10)
  x <- c(1,2,4,NA,7)
  expected_irank <- c(8, 7.5, 5, NA, 2.5)
  
  expect_error(irank_against(x,x,omega=0.5,na.rm=FALSE))
  expect_equal(irank_against(x,v,omega=0.5,na.rm=FALSE),
               expected_irank)
  expect_equal(irank_against(x,v,omega = 0.5, na.rm=TRUE),
               expected_irank[-4])
})

test_that("frank works", {
  x_1 <- c(4,4,4,3,1,10,7,7)
  expected_output <- c(0.475, 0.475, 0.475, 0.250, 0.125, 1.000, 0.800, 0.800)
  actual_output <- frank(x_1, omega=0.4, increasing = TRUE, na.rm=FALSE)
  expect_equal(actual_output, expected_output)
  
  expected_output_om0 <- c(0.375, 0.375, 0.375, 0.250, 0.125, 1.000, 0.750, 0.750)
  actual_output <- frank(x_1, omega=0, increasing = TRUE, na.rm=FALSE)
  expect_equal(actual_output, expected_output_om0)
  
  expected_output_om1 <- c(0.625, 0.625, 0.625, 0.25, 0.125, 1.0, 0.875, 0.875)
  actual_output <- frank(x_1, omega=1, increasing = TRUE, na.rm=FALSE)
  expect_equal(actual_output, expected_output_om1)
})

test_that("frank_against v argument works", {
  v <- c(1,3,4,4,4,7,7,10)
  x <- c(0,1,2,3,4,5,7,8,10,11)
  expected_irank <- c(8.5, 8, 7.5, 7, 5, 3.5, 2.5, 1.5, 1, 0.5) / 8
  expect_equal(frank_against(x, v, omega = 0.5),
               expected_irank)
})