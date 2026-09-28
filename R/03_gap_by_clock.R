# The NHANES version of Graf et al.'s Web Table 9: the Black-White difference
# in each aging measure, for the same people, in years and as Cohen's d.
#
# Graf et al. 2022 (Am J Epidemiol 191:613), Web Tables 1 and 9, HRS 2016.
# Their numbers are in data_processed/graf_2022_reference.csv, copied by hand
# from the published supplement; residuals there are the VBS-DNAm subsample.

suppressMessages({
  library(dplyr)
  library(survey)
})
source("R/utils_phenoage.R")
options(survey.lonely.psu = "adjust")

full <- read.csv("data_processed/dnam_sample.csv")
a <- full[full$analytic, ]

# Clock ages are residualized on chronological age within the analytic
# sample, as Graf et al. did. DunedinPoAm is already a rate (years of
# biological aging per calendar year) and is left as-is, as they did.
MEASURES <- data.frame(
  measure = c("blood_phenoage", "HorvathAge", "HannumAge", "PhenoAge",
              "GrimAgeMort", "GrimAge2Mort", "DunedinPoAm"),
  label = c("PhenoAge (blood chemistry)", "Horvath clock", "Hannum clock",
            "PhenoAge (epigenetic)", "GrimAge", "GrimAge2", "DunedinPoAm"),
  residualize = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)
)
for (i in seq_len(nrow(MEASURES))) {
  m <- MEASURES$measure[i]
  a[[paste0(m, "_adv")]] <- if (MEASURES$residualize[i]) age_residual(a[[m]], a$age) else a[[m]]
}

# The design is built on every participant with usable methylation data and
# then subset, so standard errors reflect NHANES's full sampling structure.
full <- left_join(full, a[, c("SEQN", grep("_adv$", names(a), value = TRUE))], by = "SEQN")
des <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTDN4YR,
                 nest = TRUE, data = full)
des_a <- subset(des, analytic)

rows <- lapply(seq_len(nrow(MEASURES)), function(i) {
  v <- paste0(MEASURES$measure[i], "_adv")
  f <- as.formula(paste0("~", v))
  means <- svyby(f, ~race, des_a, svymean)
  fit <- svyglm(as.formula(paste0(v, " ~ I(race == 'Black')")), design = des_a)
  ci <- confint(fit, df = degf(des_a))[2, ]
  cbind(
    data.frame(
      measure = MEASURES$measure[i], label = MEASURES$label[i],
      resid_white = means[means$race == "White", v],
      resid_black = means[means$race == "Black", v],
      gap = coef(fit)[[2]], gap_lower = ci[[1]], gap_upper = ci[[2]]
    ),
    design_cohens_d(des_a, v)
  )
})
out <- do.call(rbind, rows) |>
  left_join(read.csv("data_processed/graf_2022_reference.csv"), by = "measure")

write.csv(out, "output/tables/gap_by_clock.csv", row.names = FALSE)
print(out |> mutate(across(where(is.numeric), \(x) round(x, 2))) |>
        select(label, resid_white, resid_black, gap, d, d_lower, d_upper, graf_cohens_d),
      row.names = FALSE)
