#------------------------------------------------------------------------------------------
#   Project             : Replicating the Blood vs. Epigenetic PhenoAge Divergence
#   Repository          : BloodVsMethylation
#   Release Version     : 0.1.0.0
#   Author              : Iris Ivy Gauran
#   Description         : Tests: Rubin's Rules Pooling for Multiply Imputed Estimates
#------------------------------------------------------------------------------------------

source(file.path("..", "..", "R", "utils_mi.R"), chdir = TRUE)

test_that("identical imputations pool to the single-copy answer", {
  p <- pool_rubin(rep(2.5, 5), rep(0.04, 5))
  expect_equal(p$estimate, 2.5)
  expect_equal(p$se, 0.2)
  expect_equal(p$df, Inf)
  expect_equal(p$lower, 2.5 - qnorm(0.975) * 0.2, tolerance = 1e-10)
})

test_that("between-copy spread widens the standard error by (1 + 1/m)", {
  est <- c(1, 2, 3, 4, 5)
  p <- pool_rubin(est, rep(1, 5))
  expect_equal(p$estimate, 3)
  expect_equal(p$se^2, 1 + (1 + 1 / 5) * var(est))
})

test_that("mismatched inputs fail loudly", {
  expect_error(pool_rubin(1:5, 1:4))
  expect_error(pool_rubin(1, 1))
})
