source(file.path("..", "..", "R", "utils_download.R"), chdir = TRUE)

test_that("an HTML error page is not mistaken for data", {
  f <- tempfile()
  writeLines(c("<!DOCTYPE html>", strrep("x", 20000)), f)
  expect_false(looks_like_data(f))
})

test_that("a truncated file is not mistaken for data", {
  f <- tempfile()
  writeBin(as.raw(rep(65L, 100)), f)
  expect_false(looks_like_data(f))
})

test_that("a plausible data file passes", {
  f <- tempfile()
  writeBin(charToRaw(paste0("HEADER RECORD", strrep(" ", 20000))), f)
  expect_true(looks_like_data(f))
})
