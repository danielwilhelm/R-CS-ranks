context("compare marginal and simultaneous CS for ranks")
V <- diag(rep(1,5))
test_that("simultaneous CS is not wider than marginal CS", {
  for (cstype in c("two-sided", "lower", "upper")) {
    for (stepdown in c(TRUE, FALSE)) {
      set.seed(100)
      res1M <- csranks_marg(1:5, V, coverage = 0.95, cstype = cstype, stepdown = stepdown, R = 10000, indices = 1:5)
      set.seed(100)
      res1S <- csranks_simul(1:5, V, coverage = 0.95, cstype = cstype, stepdown = stepdown, R = 10000, indices = 1:5)
      set.seed(100)
      res2M <- csranks_marg(1:5, V, coverage = 0.95, cstype = cstype, stepdown = stepdown, R = 10000, indices = 1)
      set.seed(100)
      res2S <- csranks_simul(1:5, V, coverage = 0.95, cstype = cstype, stepdown = stepdown, R = 10000, indices = 1)
      set.seed(100)
      res3M <- csranks_marg(1:5, V, coverage = 0.95, cstype = cstype, stepdown = stepdown, R = 10000, indices = 1:3)
      set.seed(100)
      res3S <- csranks_simul(1:5, V, coverage = 0.95, cstype = cstype, stepdown = stepdown, R = 10000, indices = 1:3)

      expect_true(all(res1M$L >= res1S$L & res1M$U <= res1S$U))
      expect_true(all(res2M$L >= res2S$L & res2M$U <= res2S$U))
      expect_true(all(res3M$L >= res3S$L & res3M$U <= res3S$U))
    }
  }
})


test_that("CS with stepdown is not wider than single-step CS", {
  for (cstype in c("two-sided", "lower", "upper")) {
    res1T <- csranks_marg(1:5, V, coverage = 0.95, cstype = cstype, 
                          stepdown = TRUE, R = 1000, indices = 1:5, seed = 100)
    res1F <- csranks_marg(1:5, V, coverage = 0.95, cstype = cstype, 
                          stepdown = FALSE, R = 1000, indices = 1:5, seed = 100)

    expect_true(all(res1T$L >= res1F$L & res1T$U <= res1F$U))

    res1T <- csranks_simul(1:5, V, coverage = 0.95, cstype = cstype, 
                           stepdown = TRUE, R = 1000, indices = 1:5, seed = 100)
    res1F <- csranks_simul(1:5, V, coverage = 0.95, cstype = cstype, 
                           stepdown = FALSE, R = 1000, indices = 1:5, seed = 100)

    expect_true(all(res1T$L >= res1F$L & res1T$U <= res1F$U))

    res2T <- csranks_marg(1:5, V, coverage = 0.95, cstype = cstype, 
                          stepdown = TRUE, R = 1000, indices = 1, seed = 100)
    res2F <- csranks_marg(1:5, V, coverage = 0.95, cstype = cstype, 
                          stepdown = FALSE, R = 1000, indices = 1, seed = 100)

    expect_true(all(res2T$L >= res2F$L & res2T$U <= res2F$U))

    res2T <- csranks_simul(1:5, V, coverage = 0.95, cstype = cstype, 
                           stepdown = TRUE, R = 1000, indices = 1, seed = 100)
    res2F <- csranks_simul(1:5, V, coverage = 0.95, cstype = cstype, 
                           stepdown = FALSE, R = 1000, indices = 1, seed = 100)

    expect_true(all(res2T$L >= res2F$L & res2T$U <= res2F$U))

    res3T <- csranks_marg(1:5, V, coverage = 0.95, cstype = cstype, 
                          stepdown = TRUE, R = 1000, indices = 1:3, seed = 100)
    res3F <- csranks_marg(1:5, V, coverage = 0.95, cstype = cstype, 
                          stepdown = FALSE, R = 1000, indices = 1:3, seed = 100)

    expect_true(all(res3T$L >= res3F$L & res3T$U <= res3F$U))

    res3T <- csranks_simul(1:5, V, coverage = 0.95, cstype = cstype, 
                           stepdown = TRUE, R = 1000, indices = 1:3, seed = 100)
    res3F <- csranks_simul(1:5, V, coverage = 0.95, cstype = cstype, 
                           stepdown = FALSE, R = 1000, indices = 1:3, seed = 100)

    expect_true(all(res3T$L >= res3F$L & res3T$U <= res3F$U))
  }
})
