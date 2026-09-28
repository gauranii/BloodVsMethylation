#------------------------------------------------------------------------------------------
#   Project             : Replicating the Blood vs. Epigenetic PhenoAge Divergence
#   Repository          : BloodVsMethylation
#   Release Version     : 0.1.0.0
#   Author              : Iris Ivy Gauran
#   Description         : Merge Raw Files and Compute Blood-Chemistry PhenoAge
#------------------------------------------------------------------------------------------

# Join demographics, blood chemistry, blood counts, DNA-methylation clocks, and
# linked mortality into one analysis file, and compute blood PhenoAge.

suppressMessages({
  library(haven)
  library(dplyr)
})
source("R/utils_phenoage.R")

read_raw <- function(f) haven::zap_labels(haven::read_xpt(file.path("data_raw", f)))

# ---- Demographics ---------------------------------------------------------
demo <- bind_rows(read_raw("DEMO.xpt"), read_raw("DEMO_B.xpt")) |>
  transmute(
    SEQN,
    age = RIDAGEYR,
    sex = if_else(RIAGENDR == 1, "Male", "Female"),
    race = case_when(
      RIDRETH1 == 3 ~ "White",
      RIDRETH1 == 4 ~ "Black",
      RIDRETH1 %in% c(1, 2) ~ "Hispanic",
      TRUE ~ "Other"
    ),
    SDMVPSU, SDMVSTRA, WTMEC4YR
  )

# ---- Blood chemistry --------------------------------------------------------
# Alkaline phosphatase is LBXSAPSI in 1999-2000 and LBDSAPSI in 2001-2002.
chem <- bind_rows(
  read_raw("LAB18.xpt") |> transmute(SEQN, albumin = LBDSALSI, creatinine = LBDSCRSI,
                                     glucose = LBDSGLSI, alp = LBXSAPSI, iron = LBXSIR),
  read_raw("L40_B.xpt") |> transmute(SEQN, albumin = LBDSALSI, creatinine = LBDSCRSI,
                                     glucose = LBDSGLSI, alp = LBDSAPSI, iron = LBXSIR)
)

crp <- bind_rows(read_raw("LAB11.xpt"), read_raw("L11_B.xpt")) |>
  transmute(SEQN, crp_mg_dl = LBXCRP)

cbc <- bind_rows(read_raw("LAB25.xpt"), read_raw("L25_B.xpt")) |>
  transmute(SEQN, wbc = LBXWBCSI, lymph_pct = LBXLYPCT, mcv = LBXMCVSI, rdw = LBXRDW,
            hemoglobin = LBXHGB,
            neut_pct = LBXNEPCT, neut_count = LBDNENO, lymph_count = LBDLYMNO)

# ---- DNA-methylation clocks -------------------------------------------------
# A zero weight marks participants without usable methylation data; every
# clock is missing for exactly those rows.
dnam <- haven::zap_labels(haven::read_sas("data_raw/dnmepi.sas7bdat")) |>
  filter(WTDN4YR > 0)

# ---- Linked mortality (public-use, follow-up through 2019-12-31) -----------
# Fixed-width layout from the NCHS public-use LMF read-in program.
read_mort <- function(f) {
  utils::read.fwf(
    file.path("data_raw", f),
    widths = c(6, 8, 1, 1, 3, 1, 1, 21, 3, 3),
    col.names = c("SEQN", "skip1", "eligstat", "mortstat", "ucod_leading",
                  "diabetes", "hyperten", "skip2", "permth_int", "permth_exam"),
    colClasses = c("integer", "character", "integer", "integer", "character",
                   "integer", "integer", "character", "integer", "integer"),
    na.strings = c(".", "")
  ) |>
    select(SEQN, eligstat, mortstat, ucod_leading, permth_exam)
}
mort <- bind_rows(
  read_mort("NHANES_1999_2000_MORT_2019_PUBLIC.dat"),
  read_mort("NHANES_2001_2002_MORT_2019_PUBLIC.dat")
)

# ---- Join -------------------------------------------------------------------
full <- dnam |>
  inner_join(demo, by = "SEQN") |>
  left_join(chem, by = "SEQN") |>
  left_join(crp, by = "SEQN") |>
  left_join(cbc, by = "SEQN") |>
  left_join(mort, by = "SEQN") |>
  mutate(log_crp = log(crp_mg_dl))

full$blood_complete <- stats::complete.cases(full[, c(PHENOAGE_MARKERS, "age")])
full$blood_phenoage <- ifelse(full$blood_complete, phenoage_blood(full), NA_real_)
full$analytic <- full$blood_complete & full$race %in% c("White", "Black")

dir.create("data_processed", showWarnings = FALSE)
write.csv(full, "data_processed/dnam_sample.csv", row.names = FALSE)

flow <- data.frame(
  step = c("DNAm file, usable methylation (WTDN4YR > 0)",
           "  ...with all nine blood markers",
           "  ...and non-Hispanic White or Black (analytic sample)",
           "     non-Hispanic White",
           "     non-Hispanic Black"),
  n = c(nrow(full), sum(full$blood_complete), sum(full$analytic),
        sum(full$analytic & full$race == "White"),
        sum(full$analytic & full$race == "Black"))
)
dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(flow, "output/tables/sample_flow.csv", row.names = FALSE)
print(flow, row.names = FALSE)
