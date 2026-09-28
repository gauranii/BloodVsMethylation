#------------------------------------------------------------------------------------------
#   Project             : Replicating the Blood vs. Epigenetic PhenoAge Divergence
#   Repository          : BloodVsMethylation
#   Release Version     : 0.1.0.0
#   Author              : Iris Ivy Gauran
#   Description         : Red Cell Distribution Width: Iron, Hemoglobin Variants, and Mortality
#------------------------------------------------------------------------------------------

# Red cell distribution width is the largest single contributor to the blood
# PhenoAge gap (R/04_decompose_blood_gap.R). Three checks on why:
#
# 1. How much of the RDW gap survives adjustment for iron status (serum iron,
#    from the standard biochemistry panel) and other red-cell measures.
# 2. A literature-based estimate of how much inherited hemoglobin variants
#    (alpha-thalassemia -3.7 deletion, sickle cell trait, hemoglobin C trait)
#    would produce, using effect sizes and carrier frequencies from the
#    Jackson Heart Study. NHANES has no globin genotypes, so this is an
#    estimate from published numbers, not a measurement in this sample.
# 3. Whether RDW predicts death equally well in Black and White participants.

suppressMessages({
  library(dplyr)
  library(survey)
  library(survival)
})
source("R/utils_phenoage.R")
options(survey.lonely.psu = "adjust")

full <- read.csv("data_processed/dnam_sample.csv") |>
  mutate(black = as.numeric(race == "Black"),
         female = as.numeric(sex == "Female"),
         low_iron = as.numeric(iron < 60))  # ug/dL, a conventional low cutoff
des <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTDN4YR,
                 nest = TRUE, data = full)
des_iron <- subset(des, analytic & !is.na(iron) & !is.na(hemoglobin))
rdw_years <- phenoage_years_per_unit()[["rdw"]]
mcv_years <- phenoage_years_per_unit()[["mcv"]]

# ---- 1. Adjustment ----------------------------------------------------------
# The last model adjusts for MCV, which alpha-thalassemia also lowers, so it
# removes some of the inherited part along with iron; the rows are nested
# views of the same gap, not pieces to add up.
MODELS <- c(
  "Age and sex" = "rdw ~ age + female + black",
  "+ serum iron (log)" = "rdw ~ age + female + black + log(iron)",
  "+ serum iron, hemoglobin, MCV" = "rdw ~ age + female + black + log(iron) + hemoglobin + mcv"
)
adjust <- do.call(rbind, lapply(names(MODELS), function(m) {
  fit <- svyglm(as.formula(MODELS[[m]]), design = des_iron)
  ci <- confint(fit, "black", df = degf(des_iron))
  data.frame(adjustment = m, rdw_gap = coef(fit)[["black"]],
             lower = ci[1], upper = ci[2],
             phenoage_years = coef(fit)[["black"]] * rdw_years)
}))
low_iron <- svyby(~low_iron, ~race, des_iron, svymean)

# ---- 2. Inherited hemoglobin variants (literature-based) --------------------
globin <- read.csv("data_processed/globin_variant_reference.csv") |>
  mutate(
    rdw_gap_expected = (mean_copies_black - mean_copies_white) * rdw_effect_per_copy,
    mcv_gap_expected = (mean_copies_black - mean_copies_white) * coalesce(mcv_effect_per_copy, 0),
    phenoage_years_via_rdw = rdw_gap_expected * rdw_years,
    phenoage_years_via_mcv = mcv_gap_expected * mcv_years,
    phenoage_years_net = phenoage_years_via_rdw + phenoage_years_via_mcv
  )
observed <- read.csv("output/tables/blood_phenoage_decomposition.csv")
globin_summary <- data.frame(
  quantity = c("RDW gap (points)", "MCV gap (fL)", "PhenoAge years via RDW",
               "PhenoAge years via MCV", "PhenoAge years, net"),
  expected_from_variants = c(sum(globin$rdw_gap_expected), sum(globin$mcv_gap_expected),
                             sum(globin$phenoage_years_via_rdw), sum(globin$phenoage_years_via_mcv),
                             sum(globin$phenoage_years_net)),
  observed = c(observed$marker_gap_age_adjusted[observed$marker == "rdw"],
               observed$marker_gap_age_adjusted[observed$marker == "mcv"],
               observed$contribution_years[observed$marker == "rdw"],
               observed$contribution_years[observed$marker == "mcv"],
               sum(observed$contribution_years[observed$marker %in% c("rdw", "mcv")]))
)

# ---- 3. RDW and mortality by race -------------------------------------------
a <- full[full$analytic, ]
full$rdw_z <- (full$rdw - mean(a$rdw)) / sd(a$rdw)
des_m <- subset(svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTDN4YR,
                          nest = TRUE, data = full),
                analytic & eligstat == 1 & !is.na(permth_exam))
hr_by_race <- do.call(rbind, lapply(c("White", "Black"), function(g) {
  fit <- svycoxph(Surv(permth_exam, mortstat) ~ age + female + rdw_z,
                  design = subset(des_m, race == g))
  ci <- exp(confint(fit)["rdw_z", ])
  data.frame(race = g, hr_per_sd = exp(coef(fit)[["rdw_z"]]), lower = ci[[1]], upper = ci[[2]])
}))
inter <- svycoxph(Surv(permth_exam, mortstat) ~ age + female + black * rdw_z, design = des_m)
hr_by_race$interaction_p <- summary(inter)$coefficients["black:rdw_z", ncol(summary(inter)$coefficients)]

write.csv(adjust, "output/tables/rdw_adjustment.csv", row.names = FALSE)
write.csv(data.frame(race = low_iron$race, share_low_serum_iron = low_iron$low_iron),
          "output/tables/rdw_low_iron.csv", row.names = FALSE)
write.csv(globin, "output/tables/rdw_globin_variants.csv", row.names = FALSE)
write.csv(globin_summary, "output/tables/rdw_globin_summary.csv", row.names = FALSE)
write.csv(hr_by_race, "output/tables/rdw_mortality_by_race.csv", row.names = FALSE)

rnd <- function(x) mutate(x, across(where(is.numeric), \(v) round(v, 3)))
print(rnd(adjust), row.names = FALSE)
print(rnd(as.data.frame(low_iron)), row.names = FALSE)
print(rnd(globin_summary), row.names = FALSE)
print(rnd(hr_by_race), row.names = FALSE)
