V <- diag(1:10)
V_with_NAs <- diag(c(1:8, NA, 10))

expected_args <- list(x = c(1:8, 10),
                      Sigma = diag(c(1:8, 10)),
                      indices = 1:9)

test_that("`x` argument is handled correctly by process_csranks_args",{
  expect_equal(process_csranks_args(1:10, V, NA, na.rm=TRUE),
               list(x=1:10, Sigma=V, indices=1:10))
  expect_error(process_csranks_args(data.frame(x=1:10),
                                    V, NA, na.rm=TRUE))
  expect_error(process_csranks_args(paste0(1:10),
                                    V, NA, na.rm=TRUE))
})

test_that("`Sigma` argument is handled correctly by process_csranks_args",{
  expect_error(process_csranks_args(1:10, data.frame(x=1:10), NA, na.rm=TRUE))
  expect_error(process_csranks_args(1:10, 1:10 / 10, NA, na.rm=TRUE),
               "variances")
  expect_error(process_csranks_args(1:10, V[-1,], NA, na.rm=TRUE))
  expect_error(process_csranks_args(1:10, V[,-1], NA, na.rm=TRUE))
})

test_that("NAs are handled correctly by process_csranks_args", {
  expect_error(process_csranks_args(c(1:8, NA, 10),
                                    V, NA, na.rm=FALSE))
  expect_error(process_csranks_args(c(1:10),
                                    V_with_NAs, NA, na.rm=FALSE))
  expect_equal(process_csranks_args(c(1:8, NA, 10),
                                    V, NA, na.rm=TRUE),
               expected_args)
  expect_equal(process_csranks_args(1:10,
                                    V_with_NAs, NA, na.rm=TRUE),
               expected_args)
})

test_that("NAs are handled correctly at the higher level", {
  expect_error(csranks(c(1:8, NA, 10), V, na.rm=FALSE))
  expect_error(csranks(c(1:10), V_with_NAs, na.rm=FALSE))
  out_handled_NAs_internally <- csranks(c(1:8, NA, 10), V, na.rm=TRUE, seed=2023)
  out_handled_NAs_manually <- csranks(c(1:8, 10), diag(c(1:8, 10)), seed=2023)
  expect_equal(out_handled_NAs_internally, out_handled_NAs_manually)
  out_handled_NAs_internally <- csranks(1:10, V_with_NAs, na.rm=TRUE, seed=2023)
  expect_equal(out_handled_NAs_internally, out_handled_NAs_manually)
})

test_that("adjust_indices_for_NAs work correctly",{
  expect_equal(adjust_indices_for_NAs(1:10, rep(TRUE, 10)),
               1:10)
  expect_equal(adjust_indices_for_NAs(1:10, c(rep(TRUE, 4),
                                              rep(FALSE, 3), 
                                              rep(TRUE, 3))),
               1:7)
  expect_equal(adjust_indices_for_NAs(seq(2,10,2), c(TRUE, FALSE, TRUE, FALSE, FALSE,
                                              TRUE, TRUE, TRUE, FALSE, TRUE)),
               c(3,5,6))
})

test_that("`x` argument is handled correctly by process_csranks_multinom_args",{
  expect_equal(process_csranks_multinom_args(1:10, NA, na.rm=TRUE)[["x"]],
               1:10)
  expect_equal(process_csranks_multinom_args(c(NA, 2:10), NA, na.rm=TRUE)[["x"]],
               2:10)
  expect_error(process_csranks_multinom_args(data.frame(x=1:10), NA, na.rm=TRUE))
  expect_error(process_csranks_multinom_args(paste0(1:10), NA, na.rm=TRUE))
  expect_error(process_csranks_multinom_args(c(-1,2:10), NA, na.rm=TRUE))
  expect_error(process_csranks_multinom_args(c(3.14, 2:10), NA, na.rm=TRUE))
  expect_error(process_csranks_multinom_args(c(NA, 2:10), NA, na.rm=FALSE))
})

test_that("`indices` argument is handled correctly by process_csranks_multinom_args",{
  expect_equal(process_csranks_multinom_args(1:10, NA, na.rm=TRUE)[["indices"]],
               1:10)
  expect_equal(process_csranks_multinom_args(c(NA, 2:10), NA, na.rm=TRUE)[["indices"]],
               1:9)
  expect_equal(process_csranks_multinom_args(c(NA, 2:10), 1:10, na.rm=TRUE)[["indices"]],
               1:9)
  expect_equal(process_csranks_multinom_args(c(NA, 2:10), 2:10, na.rm=TRUE)[["indices"]],
               1:9)
  expect_equal(process_csranks_multinom_args(c(NA, 2,3,NA,NA,6:10), 
                                         c(1,3,5,6), na.rm=TRUE)[["indices"]],
               c(2,3))
})

test_that("`indices` argument is handled correctly by process_indices_argument",{
  expect_equal(process_indices_argument(1:10, 10),
               1:10)
  expect_equal(process_indices_argument(c(NA, 2:10), 10),
               1:10)
  expect_equal(process_indices_argument(c(-2,-5,-6), 10),
               c(1,3,4,7,8,9,10))
  expect_error(process_indices_argument(data.frame(x=1:10), 10))
  expect_error(process_indices_argument(paste0(1:10), 10))
  expect_error(process_indices_argument(c(0,2:10), 10))
  expect_error(process_indices_argument(c(1:9, 11), 10))
  expect_error(process_indices_argument(c(3.14, 2:10), 10))
  expect_error(process_indices_argument(c(0, 2:10), 10))
  expect_error(process_indices_argument(c(-1, 2:10), 10))
  expect_error(process_indices_argument(2:11, 10))
})

test_that("`x` argument is handled correctly in process_compare_args",{
  x <- c(1,3,4,4,4,7,7,10)
  omega <- 0.4
  expect_error(process_compare_args(c(x, NA), omega=omega, increasing=TRUE,
                                  na.rm = FALSE))
  expect_equal(process_compare_args(c(x, NA), omega=omega, increasing=TRUE,
                                  na.rm = TRUE),
               x)
  expect_equal(process_compare_args(x, omega=omega, increasing=FALSE,
                                  na.rm = TRUE),
               -x)
})

# simple_lmranks
test_that("`W` argument is handled correctly by process_simple_lmranks_args", {
  Y <- 1:10
  X <- 1:10
  W <- matrix(1:10, ncol = 1)
  expected_W <- cbind(rep(1, 10), W)
  expect_equal(process_simple_lmranks_args(Y, X, W, 0.4, FALSE, FALSE),
               list(Y = Y, X = X, W = expected_W))
  expect_equal(process_simple_lmranks_args(Y, X, NULL, 0.4, FALSE, FALSE),
               list(Y = Y, X = X, W = expected_W[,1,drop=FALSE]))
  expect_error(process_simple_lmranks_args(Y, X, as.vector(W), 0.4, FALSE, FALSE))
  expect_error(process_simple_lmranks_args(Y, X, "W", 0.4, FALSE, FALSE))
})

test_that("NAs are handled correctly by process_simple_lmranks_args", {
  Y <- c(1:4, NA, 6:10)
  X <- c(1:3, NA, 5:10)
  W <- matrix(c(1:5, NA, 7:10), ncol = 1)
  expected_l <- list(Y = c(1:3, 7:10),
                     X = c(1:3, 7:10),
                     W = matrix(c(rep(1, 7),
                                1:3, 7:10), ncol = 2))
  expected_l_default_W <- list(Y = c(1:3, 6:10),
                               X = c(1:3, 6:10),
                               W = matrix(rep(1, 8), ncol = 1))
  expect_error(process_simple_lmranks_args(Y, X, W, 0.4, FALSE, na.rm = FALSE))
  expect_equal(process_simple_lmranks_args(Y, X, W, 0.4, FALSE, na.rm = TRUE),
               expected_l)
  expect_equal(process_simple_lmranks_args(Y, X, NULL, 0.4, FALSE, na.rm = TRUE),
               expected_l_default_W)
})

# Low-level assert tests
test_that("assert_is_between works correctly",{
  expect_silent(assert_is_between(2:4,1:3,3:5,"x","L","U"))
  expect_error(assert_is_between(c(2,2,2),1:3,3:5,"x","L","U"))
  expect_error(assert_is_between(2:4,1:3,c(3,3,3),"x","L","U"))
})

test_that("assert_is_positive works correctly", {
  expect_silent(assert_is_positive(0.95, "x", TRUE))
  expect_error(assert_is_positive(-1, "x", TRUE))
  expect_error(assert_is_positive(0, "x", TRUE))
})

test_that("assert_is_single_probability works correctly", {
  expect_silent(assert_is_single_probability(0.95, "x"))
  expect_error(assert_is_single_probability(-1, "x"))
  expect_error(assert_is_single_probability(1.1, "x"))
})

test_that("assert_is_integer works correctly", {
  expect_silent(assert_is_integer(1:3, "x"))
  expect_silent(assert_is_integer(NA, "x", na_ok=TRUE))
  expect_error(assert_is_integer(c(1,2,3.14), "x"))
})

test_that("assert_is_numeric_vector works correctly", {
  expect_silent(assert_is_numeric_vector(1:3, "x"))
  expect_silent(assert_is_numeric_vector(NA, "x", na_ok=TRUE))
  expect_error(assert_is_numeric_vector(LETTERS[1:3], "x"))
})

test_that("assert_is_single works correctly", {
  expect_silent(assert_is_single(1, "x"))
  expect_error(assert_is_single(1:3, "x"))
})

test_that("assert_equal_length works correctly", {
  x <- 1:10
  y <- 1:10
  z <- 1:10
  expect_silent(assert_equal_length(x, y, z, names = c("x", "y", "z")))
  expect_error(assert_equal_length(x, y[1:9], z, names = c("x", "y", "z")))
})

test_that("assert_is_single_logical works correctly", {
  expect_silent(assert_is_single_logical(TRUE, "x"))
  expect_error(assert_is_single_logical(1, "x"))
  expect_error(assert_is_single_logical(NA, "x"))
})

test_that("assert_is_one_of works correctly", {
  expect_silent(assert_is_one_of("A", "x", LETTERS[1:3]))
  expect_error(assert_is_one_of(1, "x", LETTERS[1:3]))
  expect_error(assert_is_one_of("Y", "x", LETTERS[1:3]))
})

test_that("assert_has_no_NAs works correctly", {
  expect_silent(assert_has_no_NAs(1:4, "x"))
  expect_error(assert_has_no_NAs(c(NA, 1:3), "x"))
})

test_that("assert_is_vector works correctly", {
  expect_silent(assert_is_vector(1:4, "x"))
  expect_error(assert_is_vector(list(1:4), "x"))
})
