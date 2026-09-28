#------------------------------------------------------------------------------------------
#   Project             : Replicating the Blood vs. Epigenetic PhenoAge Divergence
#   Repository          : BloodVsMethylation
#   Release Version     : 0.1.0.0
#   Author              : Iris Ivy Gauran
#   Description         : Tests: PhenoAge Closed Form, Overflow, and Per-Marker Slopes
#------------------------------------------------------------------------------------------

source(file.path("..", "..", "R", "utils_phenoage.R"), chdir = TRUE)

typical <- data.frame(albumin = 42, creatinine = 80, glucose = 5.2, log_crp = log(0.2),
                      lymph_pct = 30, mcv = 90, rdw = 13, alp = 70, wbc = 6.5, age = 60)

test_that("closed form matches the published two-step formula", {
  people <- typical[rep(1, 5), ]
  people$age <- c(50, 60, 70, 80, 85)
  people$rdw <- c(12, 13, 14, 15, 16)
  expect_equal(phenoage_blood(people), phenoage_blood_two_step(people), tolerance = 1e-8)
})

test_that("closed form stays finite where the two-step form overflows", {
  sick <- typical
  sick$creatinine <- 1200; sick$rdw <- 25; sick$glucose <- 25; sick$age <- 85
  expect_false(is.finite(phenoage_blood_two_step(sick)))
  expect_true(is.finite(phenoage_blood(sick)))
})

test_that("PhenoAge is linear in each marker with the documented slope", {
  per_unit <- phenoage_years_per_unit()
  for (k in PHENOAGE_MARKERS) {
    bumped <- typical
    bumped[[k]] <- bumped[[k]] + 1
    expect_equal(phenoage_blood(bumped) - phenoage_blood(typical), per_unit[[k]],
                 tolerance = 1e-10, info = k)
  }
})

test_that("a typical 60-year-old lands near their chronological age", {
  expect_gt(phenoage_blood(typical), 45)
  expect_lt(phenoage_blood(typical), 70)
})
