#------------------------------------------------------------------------------------------
#   Project             : Replicating the Blood vs. Epigenetic PhenoAge Divergence
#   Repository          : BloodVsMethylation
#   Release Version     : 0.1.0.0
#   Author              : Iris Ivy Gauran
#   Description         : Simulate a Pure Level Shift Against Graf et al.'s Precision Checks
#------------------------------------------------------------------------------------------

# Can Graf et al.'s two precision checks detect a pure level shift?
#
# Their checks: (a) biological-age SDs are similar in Black and White
# participants, and (b) race-by-measure interaction tests on health outcomes
# look like the same tests run on chronological age. Both are about spread and
# slope. Neither looks at the mean.
#
# The test here: copy the White participants into a synthetic second group
# with the same people, the same ages, and the same deaths, then shift that
# group's markers by the Black-White differences observed in
# R/04_decompose_blood_gap.R. Biologically the two groups are identical by
# construction. Any gap that appears is measurement, and whatever the checks
# say about it, they say it about a gap that is 100% artifact.

suppressMessages({
  library(dplyr)
  library(survival)
})
source("R/utils_phenoage.R")

full <- read.csv("data_processed/dnam_sample.csv")
a <- full[full$analytic, ]
dec <- read.csv("output/tables/blood_phenoage_decomposition.csv")
shift <- setNames(dec$marker_gap_age_adjusted, dec$marker)

white <- a[a$race == "White" & a$eligstat == 1, ]

run_scenario <- function(markers) {
  real <- mutate(white, group = "Original")
  copy <- mutate(white, group = "Shifted copy")
  for (k in markers) copy[[k]] <- copy[[k]] + shift[[k]]
  sim <- bind_rows(real, copy)
  sim$adv <- age_residual(phenoage_blood(sim), sim$age)
  sim$adv_z <- as.numeric(scale(sim$adv))
  sim$shifted <- as.numeric(sim$group == "Shifted copy")
  sim$female <- as.numeric(sim$sex == "Female")

  hr_in <- function(g) {
    fit <- coxph(Surv(permth_exam, mortstat) ~ age + female + adv_z,
                 data = sim[sim$group == g, ], weights = WTDN4YR, robust = TRUE)
    exp(coef(fit)[["adv_z"]])
  }
  inter <- coxph(Surv(permth_exam, mortstat) ~ age + female + shifted * adv_z,
                 data = sim, weights = WTDN4YR, cluster = SEQN)
  s <- summary(inter)$coefficients

  data.frame(
    shifted_markers = paste(markers, collapse = " + "),
    gap_years = weighted_mean_gap(sim$adv, sim$group, sim$WTDN4YR, "Shifted copy", "Original"),
    sd_original = sd(sim$adv[sim$group == "Original"]),
    sd_shifted = sd(sim$adv[sim$group == "Shifted copy"]),
    hr_original = hr_in("Original"),
    hr_shifted = hr_in("Shifted copy"),
    # Zero to machine precision by construction, so its p-value is rounding
    # noise over rounding noise and is not reported. The group main effect
    # absorbs the shift exactly (it equals -slope x shift in SD units), which
    # is why an outcome model that adjusts for race cannot see a level shift.
    interaction_coef = s["shifted:adv_z", "coef"],
    group_main_effect = s["shifted", "coef"]
  )
}

out <- bind_rows(
  run_scenario("creatinine"),
  run_scenario(c("rdw", "mcv")),
  run_scenario(PHENOAGE_MARKERS)
)
write.csv(out, "output/tables/level_shift_simulation.csv", row.names = FALSE)
print(out |> mutate(across(where(is.numeric), \(x) round(x, 3))), row.names = FALSE)
