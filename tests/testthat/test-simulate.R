# Tests for simulate.r.
#
# simulate.r is a top-level script that reads commandArgs() itself, so it
# can't just be source()'d inside the test process -- these tests invoke it
# as a real `Rscript simulate.r ...` subprocess (see run_rscript() in
# tests/testthat/helper-fixtures.R) inside a sandboxed copy of the expected
# raw/ + data/ + code/ layout.
#
# inParams() (called from simulate.r's very first lines) always needs an
# external genome GFF file that is not part of this repo's four supplied
# CSVs, regardless of how small the population/generation counts are. Tests
# that need a full successful run skip themselves (with an informative
# message) when no data/*.gff is present, rather than failing; the
# argument-handling and missing-dependency tests below don't need one and
# always run.

library(testthat)

test_that("simulate.r starts up and prints its environment-setup message", {
  sb <- build_sandbox()
  res <- run_rscript(sb, args = character(0))
  expect_true(any(grepl("set environment", res$output)))
})

test_that("simulate.r replaces ALL cli args with defaults when fewer than 5 are supplied (simulate.r:10-11)", {
  sb <- build_sandbox()
  # simulate.r's header comment (line 5) documents one optional argument,
  # but the code (lines 10-11) only checks length(argv) < 5 and, if so,
  # replaces the WHOLE argv with hard-coded defaults -- so even though a
  # 3rd arg ("99") is supplied here, it never takes effect: the default
  # replicate/scenario-row numbers (1, 7) show up in the startup line
  # instead of "99".
  res <- run_rscript(sb, args = c("whatever.csv", "whatever-seed.csv", "99"))
  expect_true(any(grepl("environment\\s+1\\s+-\\s+7", res$output)))
})

test_that("simulate.r fails (rather than degrading gracefully) when no genome GFF is present (hidden data/*.gff dependency)", {
  sb <- build_sandbox()
  gffs <- list.files(file.path(sb, "data"), pattern = "\\.gff$")
  skip_if(length(gffs) > 0, "A real genome GFF is present; the missing-dependency path isn't reachable here.")
  res <- run_rscript(sb, args = character(0))
  expect_false(res$status == 0)
})

test_that("a tiny end-to-end simulate.r run succeeds and writes a well-formed .rda when a real genome GFF is available", {
  sb <- build_sandbox()
  gffs <- list.files(file.path(sb, "data"), pattern = "\\.gff$")
  skip_if(length(gffs) == 0, "No genome GFF found in data/ next to raw/ -- add one to exercise the full pipeline.")

  real_input <- read.csv(file.path(sb, "raw", "input.csv"), stringsAsFactors = FALSE)
  ref_val <- real_input$Value[real_input$Type == "ref genome"]
  small_input_path <- file.path(sb, "raw", "input-small.csv")
  make_small_input_csv(small_input_path, ref_val) # keeps this a fast smoke test either way

  res <- run_rscript(
    sb,
    args = c("../raw/input-small.csv", "../raw/seed.csv", "1", "../raw/scenario.csv", "7"),
    timeout = 120
  )
  expect_equal(res$status, 0)
  expect_true(any(grepl("simulation completed", res$output)))

  rda_path <- file.path(sb, "data", "tPn--1_7.rda")
  expect_true(file.exists(rda_path))
  e <- new.env()
  load(rda_path, envir = e)
  expect_true(all(c("rec.host", "rec.transposon", "rec.offspring") %in% ls(e)))
  expect_equal(nrow(e$rec.host), 3) # generation number (2) + 1, from input-small.csv
})

test_that("choosing 'normal' for the transposon-titre distribution crashes initialization (func.r:58 + simulate.r:46)", {
  # 'normal' and 'uniform' are documented as equally valid choices for
  # "transposon population size per genome distribution" (see input.csv's
  # own Note column), but only poisson/negbin behave sensibly there (see
  # the dedicated rNumVec/inParams test in test-func.R). With 'normal' and
  # sd = 0, every host's titre rounds to exactly 0, and simulate.r:46
  # (`for(i0 in 1:length(lOc))`) then iterates over c(1, 0) instead of zero
  # times, feeding tPn.act() an out-of-range NA transposon -- which errors.
  sb <- build_sandbox()
  gffs <- list.files(file.path(sb, "data"), pattern = "\\.gff$")
  skip_if(length(gffs) == 0, "No genome GFF found in data/ -- add one to exercise this end-to-end path.")

  real_input <- read.csv(file.path(sb, "raw", "input.csv"), stringsAsFactors = FALSE)
  ref_val <- real_input$Value[real_input$Type == "ref genome"]
  bad_input_path <- file.path(sb, "raw", "input-normal-dist.csv")
  make_small_input_csv(bad_input_path, ref_val)
  df <- read.csv(bad_input_path, stringsAsFactors = FALSE)
  df$Value[df$Type == "transposon population size per genome distribution"] <- "normal"
  write.csv(df, bad_input_path, row.names = FALSE)

  res <- run_rscript(
    sb,
    args = c("../raw/input-normal-dist.csv", "../raw/seed.csv", "1", "../raw/scenario.csv", "7"),
    timeout = 60
  )
  expect_false(res$status == 0)
})
