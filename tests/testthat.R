# Run from the repo root: nix develop --command Rscript tests/testthat.R
library(testthat)
test_dir("tests/testthat", stop_on_failure = TRUE)
