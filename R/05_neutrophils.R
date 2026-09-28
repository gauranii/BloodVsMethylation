# The Duffy-null route, checked as far as public data allows. NHANES has no
# Duffy genotype, so this looks for its known footprint instead: a lower
# absolute neutrophil count (ANC) and a larger share of people below the
# conventional 1.5 x 10^3/uL neutropenia threshold, among Black participants.
# The methylation-estimated neutrophil proportion (NeuPP) is a second,
# assay-independent read on the same thing.

suppressMessages({
  library(dplyr)
  library(survey)
})
options(survey.lonely.psu = "adjust")

full <- read.csv("data_processed/dnam_sample.csv") |>
  mutate(anc_below_1_5 = as.numeric(neut_count < 1.5))
des <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTDN4YR,
                 nest = TRUE, data = full)
des_a <- subset(des, analytic)

vars <- c("wbc", "neut_count", "lymph_count", "neut_pct", "lymph_pct", "NeuPP", "anc_below_1_5")
labels <- c("WBC (10^3/uL)", "Neutrophils, absolute (10^3/uL)", "Lymphocytes, absolute (10^3/uL)",
            "Neutrophils, % of WBC", "Lymphocytes, % of WBC",
            "Neutrophil proportion, methylation-estimated", "Share with ANC < 1.5")

out <- do.call(rbind, lapply(seq_along(vars), function(i) {
  f <- as.formula(paste0("~", vars[i]))
  m <- svyby(f, ~race, des_a, svymean, na.rm = TRUE)
  data.frame(measure = labels[i],
             white = m[m$race == "White", vars[i]],
             black = m[m$race == "Black", vars[i]])
}))

write.csv(out, "output/tables/neutrophils_by_race.csv", row.names = FALSE)
print(out |> mutate(across(where(is.numeric), \(x) round(x, 3))), row.names = FALSE)
