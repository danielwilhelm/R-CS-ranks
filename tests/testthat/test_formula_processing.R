### Formula Processing Tests ###

test_that("process_lmranks_formula catches illegal formulas", {
  expect_error(process_lmranks_formula("y ~ x + w"))
  expect_error(process_lmranks_formula(r(y) ~ r(x) + r(w)))
  expect_error(process_lmranks_formula(r(y) ~ r(x) * w))
  expect_error(process_lmranks_formula(r(y) ~ r(x) + r(x):w + w))
  expect_error(process_lmranks_formula(r(y) ~ r(x):w:z + w))
  expect_error(process_lmranks_formula(r(y) ~ r(x):w + z))

  expect_silent(process_lmranks_formula(r(y) ~ r(x) + w))
  expect_silent(process_lmranks_formula(r(y) ~ (r(x) + w):G))
  expect_silent(process_lmranks_formula(r(y) ~ r(x):G))
  expect_silent(process_lmranks_formula(r(y) ~ r(x):G + G))
  expect_silent(process_lmranks_formula(r(y) ~ r(x):G - 1))
})

test_that("process_lmranks_formula returns correct indices", {
  expect_equal(
    process_lmranks_formula(r(y) ~ r(x) + w)$rank_terms_indices,
    1
  )
  expect_equal(
    process_lmranks_formula(r(y) ~ w * z + r(x))$rank_terms_indices,
    3
  )
  expect_equal(
    process_lmranks_formula(r(y) ~ w + z + w:z + r(x))$rank_terms_indices,
    3
  )
  expect_equal(
    process_lmranks_formula(r(y) ~ w * z + r(x) - z)$rank_terms_indices,
    2
  )
  expect_equal(
    process_lmranks_formula(r(y) ~ w * z)$rank_terms_indices,
    integer(0)
  )

  expect_equal(
    process_lmranks_formula(r(y) ~ (r(x) + w):G - 1)$rank_terms_indices,
    1
  )
  expect_equal(
    process_lmranks_formula(r(y) ~ (r(x) + w):G)$rank_terms_indices,
    2
  )
  expect_equal(
    process_lmranks_formula(r(y) ~ (r(x) + w):G + G)$rank_terms_indices,
    2
  )
})

test_that("process_lmranks_formula returns correct ranked_response flag", {
  expect_true(process_lmranks_formula(r(y) ~ r(x) + w)$ranked_response)
  expect_false(process_lmranks_formula(y ~ r(x) + w)$ranked_response)
})

test_that("process_lmranks_formula returns corrected formula", {
  expect_equal(
    process_lmranks_formula(r(y) ~ r(x) + w:G)$formula,
    r(y) ~ r(x) + w:G
  )
  expect_equal(
    process_lmranks_formula(r(y) ~ (r(x) + w):G - 1)$formula,
    r(y) ~ (r(x) + w):G - 1
  )
  expect_equal(
    process_lmranks_formula(r(y) ~ (r(x) + w):G)$formula,
    r(y) ~ r(x):G + w:G + G - 1
  )
  expect_equal(
    process_lmranks_formula(r(y) ~ (r(x) + w):G + G - 1)$formula,
    r(y) ~ (r(x) + w):G + G - 1
  )
  expect_equal(
    process_lmranks_formula(r(y) ~ (r(x) + w):G + G)$formula,
    r(y) ~ r(x):G + w:G + G - 1
  )
})

test_that("process_lmranks_formula env to formula", {
  env <- new.env()

  actual <- process_lmranks_formula(r(y) ~ r(x) + w, env)$formula
  expect_equal(
    environment(actual),
    env
  )

  actual <- process_lmranks_formula(r(y) ~ r(x):G, env)$formula
  expect_equal(
    environment(actual),
    env
  )

  actual <- process_lmranks_formula(r(y) ~ r(x):G - 1, env)$formula
  expect_equal(
    environment(actual),
    env
  )
})

test_that("process_lmranks_formula returns correct index for simplest fits", {
  expect_equal(
    process_lmranks_formula(r(y) ~ r(x) - 1),
    list(
      rank_terms_indices = 1,
      ranked_response = TRUE,
      formula = r(y) ~ r(x) - 1
    )
  )
  expect_equal(
    process_lmranks_formula(r(y) ~ r(x)),
    list(
      rank_terms_indices = 1,
      ranked_response = TRUE,
      formula = r(y) ~ r(x)
    )
  )
})

test_that("prohibit_interactions works", {
  expect_error(prohibit_interactions(r(y) ~ r(x):w, 2))
  expect_error(prohibit_interactions(r(y) ~ r(x) + r(x):w, 2))
  expect_error(prohibit_interactions(r(y) ~ r(x):w:z, 2))
})
