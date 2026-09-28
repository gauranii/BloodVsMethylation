#------------------------------------------------------------------------------------------
#   Project             : Replicating the Blood vs. Epigenetic PhenoAge Divergence
#   Repository          : BloodVsMethylation
#   Release Version     : 0.1.0.0
#   Author              : Iris Ivy Gauran
#   Description         : Download NHANES 1999-2002 Labs, DNAm Clocks, and Linked Mortality
#------------------------------------------------------------------------------------------

# Download every raw file the analysis needs into data_raw/, caching each one.
#
# Sources (all public, no registration):
#   - NHANES 1999-2000 and 2001-2002 demographics, standard biochemistry,
#     C-reactive protein, and complete blood count files (SAS transport .xpt)
#   - The NHANES DNA-methylation epigenetic biomarker file (dnmepi.sas7bdat),
#     released July 2024, adults 50+ from the same two cycles
#   - The NCHS public-use linked mortality files, follow-up through 2019-12-31
#   - For the creatinine check (R/11_creatinine.R): DXA body composition
#     (five multiply imputed copies per participant), surplus-serum cystatin C,
#     and body measurements (height)

source("R/utils_download.R")

NHANES_BASE <- "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public"
MORT_BASE <- "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality"

RAW_FILES <- data.frame(
  file = c(
    "DEMO.xpt", "LAB18.xpt", "LAB11.xpt", "LAB25.xpt",
    "DEMO_B.xpt", "L40_B.xpt", "L11_B.xpt", "L25_B.xpt",
    "dnmepi.sas7bdat",
    "NHANES_1999_2000_MORT_2019_PUBLIC.dat",
    "NHANES_2001_2002_MORT_2019_PUBLIC.dat",
    "dxx.xpt", "dxx_b.xpt", "SSCYST_A.xpt", "SSCYST_B.xpt", "BMX.xpt", "BMX_B.xpt"
  ),
  url = c(
    file.path(NHANES_BASE, "1999/DataFiles", c("DEMO.xpt", "LAB18.xpt", "LAB11.xpt", "LAB25.xpt")),
    file.path(NHANES_BASE, "2001/DataFiles", c("DEMO_B.xpt", "L40_B.xpt", "L11_B.xpt", "L25_B.xpt")),
    "https://wwwn.cdc.gov/nchs/data/nhanes/dnam/dnmepi.sas7bdat",
    file.path(MORT_BASE, c(
      "NHANES_1999_2000_MORT_2019_PUBLIC.dat",
      "NHANES_2001_2002_MORT_2019_PUBLIC.dat"
    )),
    "https://wwwn.cdc.gov/nchs/data/nhanes/public/1999/datafiles/dxx.xpt",
    "https://wwwn.cdc.gov/nchs/data/nhanes/public/2001/datafiles/dxx_b.xpt",
    file.path(NHANES_BASE, "1999/DataFiles", c("SSCYST_A.xpt")),
    file.path(NHANES_BASE, "2001/DataFiles", c("SSCYST_B.xpt")),
    file.path(NHANES_BASE, "1999/DataFiles", c("BMX.xpt")),
    file.path(NHANES_BASE, "2001/DataFiles", c("BMX_B.xpt"))
  ),
  stringsAsFactors = FALSE
)

options(timeout = max(600, getOption("timeout")))
invisible(mapply(download_raw, RAW_FILES$file, RAW_FILES$url))
