# Run the whole pipeline in order. Raw downloads are cached in data_raw/.
scripts <- c(
  "R/01_pull_data.R",
  "R/02_build_dataset.R",
  "R/03_gap_by_clock.R",
  "R/04_decompose_blood_gap.R",
  "R/05_neutrophils.R",
  "R/06_mortality.R",
  "R/07_level_shift_simulation.R",
  "R/08_figures.R",
  "R/09_sensitivity.R"
)
for (s in scripts) {
  message("\n==== ", s, " ====")
  source(s, local = new.env())
}
writeLines(capture.output(sessionInfo()), "output/session_info.txt")
