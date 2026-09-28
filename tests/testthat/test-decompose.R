source(file.path("..", "..", "R", "utils_phenoage.R"), chdir = TRUE)

make_people <- function(n, seed) {
  set.seed(seed)
  data.frame(albumin = rnorm(n, 42, 3), creatinine = rnorm(n, 80, 15),
             glucose = rnorm(n, 5.5, 1), log_crp = rnorm(n, -1.5, 1),
             lymph_pct = rnorm(n, 30, 7), mcv = rnorm(n, 90, 5), rdw = rnorm(n, 13, 1),
             alp = rnorm(n, 70, 20), wbc = rnorm(n, 6.5, 1.5), age = runif(n, 50, 85))
}

test_that("marker contributions add up exactly to the total gap", {
  d <- make_people(400, 1)
  group <- rep(c("White", "Black"), each = 200)
  d$creatinine[group == "Black"] <- d$creatinine[group == "Black"] + 15
  d$rdw[group == "Black"] <- d$rdw[group == "Black"] + 0.7
  w <- runif(400, 0.5, 2)
  total <- weighted_mean_gap(age_residual(phenoage_blood(d), d$age), group, w)
  expect_equal(sum(decompose_gap(d, group, w)$contribution_years), total, tolerance = 1e-10)
})

test_that("identical groups decompose to zero everywhere", {
  d <- make_people(100, 2)
  d <- rbind(d, d)
  group <- rep(c("White", "Black"), each = 100)
  dec <- decompose_gap(d, group, rep(1, 200))
  expect_equal(dec$contribution_years, rep(0, 9), tolerance = 1e-10)
})

test_that("a pure shift in one marker is attributed to that marker alone", {
  d <- make_people(100, 3)
  shifted <- d
  shifted$creatinine <- shifted$creatinine + 20
  both <- rbind(d, shifted)
  group <- rep(c("White", "Black"), each = 100)
  dec <- decompose_gap(both, group, rep(1, 200))
  expect_equal(dec$contribution_years[dec$marker == "creatinine"],
               20 * phenoage_years_per_unit()[["creatinine"]], tolerance = 1e-8)
  expect_equal(dec$contribution_years[dec$marker != "creatinine"], rep(0, 8), tolerance = 1e-8)
})
