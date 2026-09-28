#------------------------------------------------------------------------------------------
#   Project             : Replicating the Blood vs. Epigenetic PhenoAge Divergence
#   Repository          : BloodVsMethylation
#   Release Version     : 0.1.0.0
#   Author              : Iris Ivy Gauran
#   Description         : Shared Functions: PhenoAge, Age Residuals, Gap Decomposition, Cohen's d
#------------------------------------------------------------------------------------------

# Blood-chemistry Phenotypic Age and the arithmetic this repo builds on it.
#
# Coefficients are Levine et al. 2018 (Aging 10:573), Table 1, as used in
# Liu et al. 2018 (PLOS Medicine 15:e1002718). Units are the ones that paper
# specifies: albumin g/L, creatinine umol/L, glucose mmol/L, CRP mg/dL (log),
# lymphocyte %, MCV fL, RDW %, alkaline phosphatase U/L, WBC 1000 cells/uL.

PHENOAGE_COEF <- c(
  albumin    = -0.0336,
  creatinine =  0.0095,
  glucose    =  0.1953,
  log_crp    =  0.0954,
  lymph_pct  = -0.0120,
  mcv        =  0.0268,
  rdw        =  0.3306,
  alp        =  0.00188,
  wbc        =  0.0554,
  age        =  0.0804
)
PHENOAGE_INTERCEPT <- -19.907
GOMPERTZ_GAMMA <- 0.0076927

PHENOAGE_MARKERS <- setdiff(names(PHENOAGE_COEF), "age")

phenoage_xb <- function(d) {
  xb <- PHENOAGE_INTERCEPT
  for (k in names(PHENOAGE_COEF)) xb <- xb + PHENOAGE_COEF[[k]] * d[[k]]
  xb
}

# The published two-step form: 10-year mortality score, then the age at which
# that score is the population average. Kept for the test that checks it
# against phenoage_blood(); it returns Inf once the mortality score rounds to
# exactly 1 in floating point, which happens for the highest-risk participants.
phenoage_blood_two_step <- function(d) {
  xb <- phenoage_xb(d)
  mort <- 1 - exp(-exp(xb) * (exp(120 * GOMPERTZ_GAMMA) - 1) / GOMPERTZ_GAMMA)
  141.50225 + log(-0.00553 * log(1 - mort)) / 0.090165
}

# Substituting the mortality score into the age formula cancels both logs,
# so PhenoAge is exactly linear in xb: one unit of xb is 1/0.090165 = 11.09
# years. This closed form is what the repo uses (it never overflows), and the
# linearity is what makes the per-marker decomposition below exact.
YEARS_PER_XB <- 1 / 0.090165
phenoage_blood <- function(d) {
  141.50225 + (log(0.00553) + log((exp(120 * GOMPERTZ_GAMMA) - 1) / GOMPERTZ_GAMMA) +
                 phenoage_xb(d)) * YEARS_PER_XB
}
phenoage_years_per_unit <- function() PHENOAGE_COEF[PHENOAGE_MARKERS] * YEARS_PER_XB

# Graf et al. define "age advancement" as the residual from a regression of
# biological age on chronological age, fit once in the whole analysis sample.
age_residual <- function(y, age) {
  unname(stats::residuals(stats::lm(y ~ age, na.action = stats::na.exclude)))
}

weighted_mean_gap <- function(x, group, w, focal = "Black", ref = "White") {
  stats::weighted.mean(x[group == focal], w[group == focal]) -
    stats::weighted.mean(x[group == ref], w[group == ref])
}

# Black-White gap in blood PhenoAge advancement, split by marker. Because
# PhenoAge is linear in the markers and residualizing on age is a linear
# operator, advancement = sum_k years_per_unit_k * resid(marker_k), so the
# group gap is the sum of each marker's gap times its weight, exactly.
decompose_gap <- function(d, group, w, focal = "Black", ref = "White") {
  per_unit <- phenoage_years_per_unit()
  rows <- lapply(PHENOAGE_MARKERS, function(k) {
    r <- age_residual(d[[k]], d$age)
    marker_gap <- weighted_mean_gap(r, group, w, focal, ref)
    data.frame(
      marker = k,
      years_per_unit = unname(per_unit[[k]]),
      marker_gap_age_adjusted = marker_gap,
      contribution_years = unname(per_unit[[k]]) * marker_gap
    )
  })
  do.call(rbind, rows)
}

# Standardized mean difference with a design-based confidence interval: the
# outcome is scaled by its weighted SD in the analysis domain, then regressed
# on a focal-group indicator.
design_cohens_d <- function(design, var, group_var = "race", focal = "Black") {
  vars <- design$variables
  sd_w <- sqrt(survey::svyvar(stats::as.formula(paste0("~", var)), design, na.rm = TRUE)[1])
  design <- stats::update(design,
    .z = vars[[var]] / sd_w,
    .focal = as.numeric(vars[[group_var]] == focal)
  )
  fit <- survey::svyglm(.z ~ .focal, design = design)
  est <- stats::coef(fit)[[".focal"]]
  ci <- suppressWarnings(stats::confint(fit, ".focal", df = survey::degf(design)))
  data.frame(d = est, d_lower = ci[1], d_upper = ci[2],
             p = summary(fit)$coefficients[".focal", 4])
}
