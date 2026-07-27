# Root-level test runner for func.r / simulate.r.
#
# Run from the repository root with:
#   Rscript tests/testthat.R
#
# (Mirrors the standard R package convention of tests/testthat.R +
# tests/testthat/, adapted for a repo that is a loose collection of
# scripts rather than an installable package -- there is no DESCRIPTION
# file, so this calls test_dir() directly instead of test_check().)

library(testthat)

test_dir("tests/testthat", stop_on_failure = TRUE)
