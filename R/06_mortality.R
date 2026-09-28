#------------------------------------------------------------------------------------------
#   Project             : Replicating the Blood vs. Epigenetic PhenoAge Divergence
#   Repository          : BloodVsMethylation
#   Release Version     : 0.1.0.0
#   Author              : Iris Ivy Gauran
#   Description         : Mortality Prediction by Race, with Suspect Markers Neutralized
#------------------------------------------------------------------------------------------

# Two mortality questions, both with design-weighted Cox models on deaths
# through 2019-12-31 (months from the exam):
#
# 1. Does each aging measure predict death equally well in Black and White
#    participants? This is Liu et al.'s subgroup claim, and the check Graf et
#    al. used to argue measurement precision is similar across groups.
# 2. If a suspect marker is held at the sample mean (its contribution to
#    everyone's PhenoAge removed), how much of the Black-White gap goes with
#    it, and how much mortality prediction does the score lose?

suppressMessages({
  library(dplyr)
  library(survey)
  library(survival)
})
source("R/utils_phenoage.R")
options(survey.lonely.psu = "adjust")

full <- read.csv("data_processed/dnam_sample.csv")
a <- full[full$analytic, ]

# ---- Variants of blood PhenoAge with marker sets neutralized ---------------
NEUTRALIZE <- list(
  "Published (all nine markers)" = character(0),
  "Creatinine held at mean" = "creatinine",
  "WBC + lymphocyte % held at mean" = c("wbc", "lymph_pct"),
  "RDW + MCV held at mean" = c("rdw", "mcv"),
  "All five held at mean" = c("creatinine", "wbc", "lymph_pct", "rdw", "mcv")
)
neutralized_phenoage <- function(d, markers, w) {
  for (k in markers) d[[k]] <- stats::weighted.mean(d[[k]], w)
  phenoage_blood(d)
}

variant_cols <- character(0)
for (i in seq_along(NEUTRALIZE)) {
  col <- paste0("variant_", i, "_adv")
  a[[col]] <- age_residual(neutralized_phenoage(a, NEUTRALIZE[[i]], a$WTDN4YR), a$age)
  variant_cols <- c(variant_cols, col)
}

CLOCKS <- c(epigenetic_phenoage = "PhenoAge", grimage = "GrimAgeMort",
            grimage2 = "GrimAge2Mort", horvath = "HorvathAge", hannum = "HannumAge")
for (k in names(CLOCKS)) a[[paste0(k, "_adv")]] <- age_residual(a[[CLOCKS[[k]]]], a$age)
a$dunedinpoam_adv <- a$DunedinPoAm

adv_cols <- c(variant_cols, paste0(names(CLOCKS), "_adv"), "dunedinpoam_adv")
adv_labels <- c(paste("Blood PhenoAge:", names(NEUTRALIZE)),
                "Epigenetic PhenoAge", "GrimAge", "GrimAge2", "Horvath", "Hannum", "DunedinPoAm")

# Standardize each measure within the analytic sample so hazard ratios are
# per SD and comparable across measures.
for (v in adv_cols) a[[paste0(v, "_z")]] <- as.numeric(scale(a[[v]]))

full <- left_join(full, a[, c("SEQN", adv_cols, paste0(adv_cols, "_z"))], by = "SEQN") |>
  mutate(black = as.numeric(race == "Black"), female = as.numeric(sex == "Female"))
des <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTDN4YR,
                 nest = TRUE, data = full)
des_a <- subset(des, analytic & eligstat == 1 & !is.na(permth_exam))

hr_row <- function(fit, term) {
  b <- coef(fit)[[term]]
  se <- sqrt(vcov(fit)[term, term])
  c(hr = exp(b), lower = exp(b - 1.96 * se), upper = exp(b + 1.96 * se),
    p = 2 * pnorm(-abs(b / se)))
}

rows <- lapply(seq_along(adv_cols), function(i) {
  z <- paste0(adv_cols[i], "_z")
  base <- "Surv(permth_exam, mortstat) ~ age + female"
  overall <- svycoxph(as.formula(paste(base, "+ black +", z)), design = des_a)
  white <- svycoxph(as.formula(paste(base, "+", z)), design = subset(des_a, race == "White"))
  black <- svycoxph(as.formula(paste(base, "+", z)), design = subset(des_a, race == "Black"))
  inter <- svycoxph(as.formula(paste(base, "+ black *", z)), design = des_a)
  gap <- weighted_mean_gap(a[[adv_cols[i]]], a$race, a$WTDN4YR)
  data.frame(
    measure = adv_labels[i],
    gap_years_or_units = gap,
    t(setNames(hr_row(overall, z), paste0("hr_all_", c("est", "lower", "upper", "p")))),
    t(setNames(hr_row(white, z), paste0("hr_white_", c("est", "lower", "upper", "p")))),
    t(setNames(hr_row(black, z), paste0("hr_black_", c("est", "lower", "upper", "p")))),
    interaction_p = hr_row(inter, paste0("black:", z))[["p"]],
    check.names = FALSE
  )
})
out <- do.call(rbind, rows)
write.csv(out, "output/tables/mortality_by_measure.csv", row.names = FALSE)

cat("deaths:", sum(a$mortstat == 1 & a$eligstat == 1), "of", sum(a$eligstat == 1), "\n")
print(out |> mutate(across(where(is.numeric), \(x) round(x, 3))) |>
        select(measure, gap_years_or_units, hr_all_est, hr_white_est, hr_black_est, interaction_p),
      row.names = FALSE)
