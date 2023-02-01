I <- matrix(c(FALSE, FALSE, FALSE, FALSE, 
              FALSE, FALSE, TRUE, TRUE, 
              FALSE, TRUE, FALSE, FALSE, 
              FALSE, FALSE, TRUE, FALSE
              ), byrow = TRUE, nrow=4)

expected_needed_variables <- c(FALSE, TRUE, TRUE, TRUE)
expected_requested_diffrences <- matrix(c(2,1,
                                           1,2,
                                           3,2,
                                           1,3),
                                         byrow = TRUE, ncol=2)
expected_reduced_I <- list(needed_variables = expected_needed_variables,
                           requested_diffrences = expected_requested_diffrences)
test_that("reduce_I works",{
  actual_reduced_I <- reduce_I(I)
  expect_equal(actual_reduced_I, expected_reduced_I)
})

Z <- matrix(c(1:3, seq(1,10,4)), byrow=TRUE, nrow=2)
scales <- c(1,1,2,3)
expected_Z_diff <- matrix(c((Z[,2] - Z[,1]) / scales[1],
                            (Z[,1] - Z[,2]) / scales[2],
                            (Z[,3] - Z[,2]) / scales[3],
                            (Z[,1] - Z[,3]) / scales[4]),
                          nrow = 2)
test_that("calculate_scaled_diffs works", {
  actual_Z_diff <- calculate_scaled_differences_in_samples(
    Z, expected_requested_diffrences,scales)
  expect_equal(actual_Z_diff, expected_Z_diff)
})