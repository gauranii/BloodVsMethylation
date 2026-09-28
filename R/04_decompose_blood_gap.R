# Split the Black-White gap in blood PhenoAge advancement into the part each
# of the nine markers contributes, in years. Exact, not approximate: see
# decompose_gap() in R/utils_phenoage.R for why.

suppressMessages(library(dplyr))
source("R/utils_phenoage.R")

full <- read.csv("data_processed/dnam_sample.csv")
a <- full[full$analytic, ]

MARKER_LABELS <- c(
  albumin = "Albumin", creatinine = "Creatinine", glucose = "Glucose",
  log_crp = "C-reactive protein (log)", lymph_pct = "Lymphocyte %",
  mcv = "Mean cell volume", rdw = "Red cell distribution width",
  alp = "Alkaline phosphatase", wbc = "White blood cell count"
)

dec <- decompose_gap(a, a$race, a$WTDN4YR) |>
  mutate(label = MARKER_LABELS[marker]) |>
  arrange(desc(contribution_years))

total_gap <- weighted_mean_gap(age_residual(a$blood_phenoage, a$age), a$race, a$WTDN4YR)
stopifnot(abs(sum(dec$contribution_years) - total_gap) < 1e-8)

# The two groups of markers the race-bias argument is about: creatinine
# (muscle mass) and the white-cell pair (Duffy-null lowers neutrophils,
# which lowers WBC and raises lymphocyte %).
summary_rows <- data.frame(
  component = c("Total gap", "Creatinine", "White-cell markers (WBC + lymphocyte %)",
                "All other markers"),
  years = c(total_gap,
            sum(dec$contribution_years[dec$marker == "creatinine"]),
            sum(dec$contribution_years[dec$marker %in% c("wbc", "lymph_pct")]),
            sum(dec$contribution_years[!dec$marker %in% c("creatinine", "wbc", "lymph_pct")]))
)

write.csv(dec, "output/tables/blood_phenoage_decomposition.csv", row.names = FALSE)
write.csv(summary_rows, "output/tables/blood_phenoage_decomposition_summary.csv", row.names = FALSE)
print(dec |> mutate(across(where(is.numeric), \(x) round(x, 3))) |>
        select(label, years_per_unit, marker_gap_age_adjusted, contribution_years), row.names = FALSE)
print(summary_rows |> mutate(years = round(years, 2)), row.names = FALSE)
