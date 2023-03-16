
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

test_that("irank without na.rm handles NAs", {
  x <- c(1,3,4,4,4,7,7,10, NA, NA, NA)
  expect_error(irank(x))
})

test_that("irank with na.rm handles NAs", {
  x <- c(1,3,4,4,4,7,7,10, NA, NA, NA)
  expected_irank_om0 <- c(8,7,4,4,4,2,2,1)
  expected_irank_om1 <- c(8,7,6,6,6,3,3,1)
  expected_irank_om0.5 <- c(8,7,5,5,5,2.5,2.5,1)
  expect_equal(irank(x, omega = 0, na.rm=TRUE),
               expected_irank_om0)
  expect_equal(irank(x, omega = 1, na.rm=TRUE),
               expected_irank_om1)
  expect_equal(irank(x, omega = 0.5, na.rm=TRUE),
               expected_irank_om0.5)
})

test_that("create_bins works", {
  x <- 1:20
  bins <- 5
  expected_bins <- factor(rep(c("1", "2", "3", "4", "5"), each = 4))
  
  expect_equal(createbins(x, bins),
               expected_bins)
})

test_that("get_double_from_single_indices works", {
  full_double_indices <- matrix(c(
    rep(1:4, times = 4),
    rep(1:4, each = 4)
  ),
  ncol = 2
  )
  
  expect_equal(
    get_double_from_single_indices(1:16, 4),
    full_double_indices
  )
  expect_equal(
    get_double_from_single_indices(as.vector(matrix(1:16, ncol = 4)), 4),
    full_double_indices
  )
})
