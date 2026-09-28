# Does the headline (blood PhenoAge gap, epigenetic PhenoAge gap, and the
# decomposition's top contributors) survive three choices a reader could
# reasonably make differently?
#
#   1. Drop the 7 participants whose PhenoAge the published two-step formula
#      returns as infinite (an implementation that uses it, as the BioAge
#      package may, would lose them).
#   2. Winsorize every blood marker at its 1st and 99th percentiles.
#   3. Ignore the survey weights.

suppressMessages(library(dplyr))
source("R/utils_phenoage.R")

full <- read.csv("data_processed/dnam_sample.csv")
a <- full[full$analytic, ]

winsorize <- function(x) {
  q <- stats::quantile(x, c(0.01, 0.99))
  pmin(pmax(x, q[[1]]), q[[2]])
}

summarize_scenario <- function(name, d, w) {
  d$blood_phenoage <- phenoage_blood(d)
  dec <- decompose_gap(d, d$race, w)
  top <- dec[order(-abs(dec$contribution_years)), ][1:3, ]
  data.frame(
    scenario = name, n = nrow(d),
    blood_gap = weighted_mean_gap(age_residual(d$blood_phenoage, d$age), d$race, w),
    epigenetic_gap = weighted_mean_gap(age_residual(d$PhenoAge, d$age), d$race, w),
    rdw_years = dec$contribution_years[dec$marker == "rdw"],
    creatinine_years = dec$contribution_years[dec$marker == "creatinine"],
    white_cell_years = sum(dec$contribution_years[dec$marker %in% c("wbc", "lymph_pct")]),
    top_three = paste(top$marker, collapse = ", ")
  )
}

finite <- a[is.finite(phenoage_blood_two_step(a)), ]
wins <- a
for (k in PHENOAGE_MARKERS) wins[[k]] <- winsorize(wins[[k]])

out <- bind_rows(
  summarize_scenario("Main analysis", a, a$WTDN4YR),
  summarize_scenario("Drop two-step infinite values", finite, finite$WTDN4YR),
  summarize_scenario("Winsorize markers at 1%/99%", wins, wins$WTDN4YR),
  summarize_scenario("Unweighted", a, rep(1, nrow(a)))
)
write.csv(out, "output/tables/sensitivity.csv", row.names = FALSE)
print(out |> mutate(across(where(is.numeric), \(x) round(x, 2))), row.names = FALSE)
