context("marginal CS for ranks")
V <- diag(rep(1, 10))
test_that("return value is of correct class and size", {
  for (cstype in c("two-sided", "lower", "upper")) {
    for (stepdown in c(TRUE, FALSE)) {
      res1 <- csranks(1:10, V, coverage = 0.95, cstype = cstype, stepdown = stepdown, simul = FALSE, R = 100, indices = NA)
      res2 <- csranks(1:10, V, coverage = 0.95, cstype = cstype, stepdown = stepdown, simul = FALSE, R = 100, indices = 1)
      res3 <- csranks(1:10, V, coverage = 0.95, cstype = cstype, stepdown = stepdown, simul = FALSE, R = 100, indices = 1:3)

      expect_is(res1$L, "integer")
      expect_is(res1$U, "integer")
      expect_is(res1$rank, "numeric")
      expect_is(res2$L, "integer")
      expect_is(res2$U, "integer")
      expect_is(res2$rank, "numeric")
      expect_is(res3$L, "integer")
      expect_is(res3$U, "integer")
      expect_is(res3$rank, "numeric")
      
      expect_equal(length(res1$L), 10)
      expect_equal(length(res1$U), 10)
      expect_equal(length(res1$rank), 10)
      expect_equal(length(res2$L), 1)
      expect_equal(length(res2$U), 1)
      expect_equal(length(res2$rank), 1)
      expect_equal(length(res3$L), 3)
      expect_equal(length(res3$U), 3)
      expect_equal(length(res3$rank), 3)
      
      expect_false(any(is.na(res1$L)))
      expect_false(any(is.na(res1$U)))
      expect_false(any(is.na(res1$rank)))
      expect_false(any(is.na(res2$L)))
      expect_false(any(is.na(res2$U)))
      expect_false(any(is.na(res2$rank)))
      expect_false(any(is.na(res3$L)))
      expect_false(any(is.na(res3$U)))
      expect_false(any(is.na(res3$rank)))
    }
  }
})

test_that("lower and upper bounds are in the correct range of values", {
  for (cstype in c("two-sided", "lower", "upper")) {
    for (stepdown in c(TRUE, FALSE)) {
      res1 <- csranks(1:10, V, coverage = 0.95, cstype = cstype, stepdown = stepdown, simul = FALSE, R = 100, indices = NA)
      res2 <- csranks(1:10, V, coverage = 0.95, cstype = cstype, stepdown = stepdown, simul = FALSE, R = 100, indices = 1)
      res3 <- csranks(1:10, V, coverage = 0.95, cstype = cstype, stepdown = stepdown, simul = FALSE, R = 100, indices = 1:3)

      expect_true(all(res1$L <= res1$rank & res1$rank <= res1$U))
      expect_true(all(res2$L <= res2$rank & res2$rank <= res2$U))
      expect_true(all(res3$L <= res3$rank & res3$rank <= res3$U))

      expect_true(all(res1$L <= 10 & res1$U <= 10 & res1$L >= 1 & res1$U >= 1))
      expect_true(all(res2$L <= 10 & res2$U <= 10 & res2$L >= 1 & res2$U >= 1))
      expect_true(all(res3$L <= 10 & res3$U <= 10 & res3$L >= 1 & res3$U >= 1))

      if (cstype == "lower") {
        expect_equal(res1$U, rep(10, 10))
        expect_equal(res2$U, 10)
        expect_equal(res3$U, rep(10, 3))
      }
      if (cstype == "upper") {
        expect_equal(res1$L, rep(1, 10))
        expect_equal(res2$L, 1)
        expect_equal(res3$L, rep(1, 3))
      }
    }
  }
})
