path <- file.path("..", "..", "data_processed", "dnam_sample.csv")
skip_if_not(file.exists(path), "run R/02_build_dataset.R first")
d <- read.csv(path)

test_that("one row per participant, all with usable methylation data", {
  expect_equal(anyDuplicated(d$SEQN), 0)
  expect_true(all(d$WTDN4YR > 0))
  expect_false(anyNA(d$PhenoAge))
})

test_that("ages sit inside the DNAm eligibility window and NHANES top-code", {
  expect_true(all(d$age >= 50 & d$age <= 85))
})

test_that("the analytic sample is complete and two-group", {
  a <- d[d$analytic, ]
  expect_setequal(unique(a$race), c("White", "Black"))
  expect_false(anyNA(a$blood_phenoage))
  expect_true(all(is.finite(a$blood_phenoage)))
  expect_true(all(a$eligstat == 1))
})

test_that("markers are in the units Levine's coefficients expect", {
  a <- d[d$analytic, ]
  expect_true(stats::median(a$albumin) > 30 && stats::median(a$albumin) < 50)       # g/L
  expect_true(stats::median(a$creatinine) > 50 && stats::median(a$creatinine) < 130) # umol/L
  expect_true(stats::median(a$glucose) > 4 && stats::median(a$glucose) < 8)          # mmol/L
})
