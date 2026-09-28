#------------------------------------------------------------------------------------------
#   Project             : Replicating the Blood vs. Epigenetic PhenoAge Divergence
#   Repository          : BloodVsMethylation
#   Release Version     : 0.1.0.0
#   Author              : Iris Ivy Gauran
#   Description         : Creatinine: Muscle Mass (DXA) vs. Kidney Function (Cystatin C)
#------------------------------------------------------------------------------------------

# Creatinine adds 1.7 years to the blood PhenoAge gap in the main analysis
# (R/04_decompose_blood_gap.R), but only about 1.1 once extreme values are
# dropped or trimmed (R/09_sensitivity.R). Two questions:
#
# 1. How much of the Black-White creatinine difference is muscle mass and how
#    much is kidney function? Serum creatinine comes from muscle and is
#    cleared by the kidneys, so it rises with either more muscle or worse
#    kidneys. DXA appendicular lean mass (arms + legs) measures the first;
#    cystatin C, a kidney marker that does not depend on muscle, measures the
#    second. The gap is estimated in one common sample under each adjustment.
# 2. Who are the extreme values, and is their creatinine kidney disease?

suppressMessages({
  library(haven)
  library(dplyr)
  library(survey)
})
source("R/utils_phenoage.R")
source("R/utils_mi.R")
options(survey.lonely.psu = "adjust")

read_raw <- function(f) haven::zap_labels(haven::read_xpt(file.path("data_raw", f)))

# Appendicular lean mass, kg, one row per participant per imputed copy.
dxa <- bind_rows(read_raw("dxx.xpt"), read_raw("dxx_b.xpt")) |>
  transmute(SEQN, imputation = `_MULT_`,
            alm_kg = (DXDLALE + DXDRALE + DXDLLLE + DXDRLLE) / 1000)
cystatin <- bind_rows(read_raw("SSCYST_A.xpt"), read_raw("SSCYST_B.xpt")) |>
  transmute(SEQN, cystatin_c = SSCYPC)

full <- read.csv("data_processed/dnam_sample.csv") |>
  left_join(cystatin, by = "SEQN") |>
  mutate(black = as.numeric(race == "Black"), female = as.numeric(sex == "Female"))

# Common sample: analytic, with cystatin C, and with lean mass in every copy.
has_lean <- dxa |> group_by(SEQN) |> summarise(ok = all(!is.na(alm_kg)), .groups = "drop") |>
  filter(ok) |> pull(SEQN)
full$common <- full$analytic & !is.na(full$cystatin_c) & full$SEQN %in% has_lean

creat_years <- phenoage_years_per_unit()[["creatinine"]]

MODELS <- c(
  "Age and sex" = "creatinine ~ age + female + black",
  "+ lean mass" = "creatinine ~ age + female + black + alm_kg",
  "+ cystatin C" = "creatinine ~ age + female + black + log(cystatin_c)",
  "+ lean mass + cystatin C" = "creatinine ~ age + female + black + alm_kg + log(cystatin_c)"
)

# Fit one model on one imputed copy of the lean-mass data.
fit_copy <- function(formula, m, data) {
  d <- left_join(data, filter(dxa, imputation == m), by = "SEQN")
  des <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTDN4YR, nest = TRUE, data = d)
  fit <- svyglm(as.formula(formula), design = subset(des, common))
  c(est = coef(fit)[["black"]], var = vcov(fit)["black", "black"])
}

gap_rows <- function(data, sample_label) {
  do.call(rbind, lapply(names(MODELS), function(name) {
    per_copy <- sapply(1:5, function(m) fit_copy(MODELS[[name]], m, data))
    p <- pool_rubin(per_copy["est", ], per_copy["var", ])
    data.frame(sample = sample_label, adjustment = name, n = sum(data$common),
               creatinine_gap = p$estimate, lower = p$lower, upper = p$upper,
               phenoage_years = p$estimate * creat_years)
  }))
}

# ---- 1. Muscle vs. kidney ---------------------------------------------------
# Run once with everyone, and once without the participants whose creatinine
# is above 176 umol/L (2 mg/dL), the range where kidney disease dominates.
EXTREME_UMOL <- 176
no_extreme <- mutate(full, common = common & creatinine <= EXTREME_UMOL)
gaps <- bind_rows(
  gap_rows(full, "Everyone"),
  gap_rows(no_extreme, "Creatinine <= 176 umol/L")
)

# ---- 2. The extreme values ----------------------------------------------------
a <- full[full$analytic, ]
extreme <- a |>
  mutate(extreme = creatinine > EXTREME_UMOL) |>
  group_by(race, extreme) |>
  summarise(n = n(),
            median_creatinine = median(creatinine),
            median_cystatin_c = median(cystatin_c, na.rm = TRUE),
            n_with_cystatin = sum(!is.na(cystatin_c)),
            .groups = "drop")

# What the extreme values do to creatinine's share of the PhenoAge gap. The
# decomposition's residualization is rerun on the reduced sample, so this is
# the same quantity R/04 reports, just without them.
dec_all <- decompose_gap(a, a$race, a$WTDN4YR)
b <- a[a$creatinine <= EXTREME_UMOL, ]
dec_trim <- decompose_gap(b, b$race, b$WTDN4YR)
contribution <- data.frame(
  sample = c("Everyone", "Creatinine <= 176 umol/L"),
  n = c(nrow(a), nrow(b)),
  creatinine_years = c(dec_all$contribution_years[dec_all$marker == "creatinine"],
                       dec_trim$contribution_years[dec_trim$marker == "creatinine"]),
  total_gap_years = c(sum(dec_all$contribution_years), sum(dec_trim$contribution_years))
)

write.csv(gaps, "output/tables/creatinine_muscle_kidney.csv", row.names = FALSE)
write.csv(extreme, "output/tables/creatinine_extreme_values.csv", row.names = FALSE)
write.csv(contribution, "output/tables/creatinine_contribution_trimmed.csv", row.names = FALSE)

rnd <- function(x) mutate(x, across(where(is.numeric), \(v) round(v, 2)))
print(rnd(gaps), row.names = FALSE)
print(rnd(extreme), row.names = FALSE)
print(rnd(contribution), row.names = FALSE)
