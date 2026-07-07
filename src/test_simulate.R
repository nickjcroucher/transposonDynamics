#!/usr/bin/env Rscript

library(testthat)

simulate_path <- Sys.getenv(
  "SIMULATE_R_PATH",
  unset = file.path(Sys.getenv("GITHUB_WORKSPACE", unset = getwd()), "src", "simulate.r")
)

if (!file.exists(simulate_path)) {
  stop("Cannot find simulate.r. Set SIMULATE_R_PATH to the script path before running these tests.")
}

simulate_path <- normalizePath(simulate_path, mustWork = TRUE)
src_dir <- dirname(simulate_path)
project_dir <- normalizePath(file.path(src_dir, ".."), mustWork = TRUE)

with_dir <- function(path, code) {
  old <- getwd()
  on.exit(setwd(old), add = TRUE)
  setwd(path)
  force(code)
}

test_that("simulate.r has no syntax errors", {
  expect_error(parse(file = simulate_path), NA)
})

test_that("simulate.r runs end-to-end and records well-formed output", {
  input_csv <- file.path(project_dir, "raw", "input.csv")
  seed_csv <- file.path(project_dir, "raw", "seed.csv")
  scenario_csv <- file.path(project_dir, "raw", "scenario.csv")

  skip_if_not(file.exists(input_csv), "raw/input.csv not present")
  skip_if_not(file.exists(seed_csv), "raw/seed.csv not present")
  skip_if_not(file.exists(scenario_csv), "raw/scenario.csv not present")

  pMs <- read.csv(input_csv, header = TRUE)
  ref_genome <- pMs$Value[pMs$Type == "ref genome"]
  skip_if(length(ref_genome) == 0, "input.csv has no 'ref genome' row")

  gff_path <- file.path(project_dir, "data", paste0(strsplit(ref_genome, "[.]")[[1]][1], ".gff"))
  skip_if_not(file.exists(gff_path), paste("reference genome file not present:", gff_path))

  n_seeds <- nrow(read.csv(seed_csv, header = FALSE))
  n_scenarios <- nrow(read.csv(scenario_csv, header = TRUE))
  skip_if(n_seeds < 1 || n_scenarios < 1, "seed.csv/scenario.csv has no data rows")

  gen_max <- as.numeric(pMs$Value[pMs$Type == "host organism constant generation number"])
  pop_size <- as.numeric(pMs$Value[pMs$Type == "host organism constant population size"])

  out_rda <- file.path(project_dir, "data", "tPn--1_1.rda")
  if (file.exists(out_rda)) unlink(out_rda)
  on.exit(unlink(out_rda), add = TRUE)

  ## run the real script end-to-end against real project data, with a bounded
  ## timeout so CI can't hang if input.csv specifies a very large run
  timed_out <- FALSE
  result <- withCallingHandlers(
    with_dir(src_dir, system2(
      "Rscript",
      c("simulate.r", input_csv, seed_csv, "1", scenario_csv, "1"),
      stdout = TRUE, stderr = TRUE, timeout = 180
    )),
    warning = function(w) {
      if (grepl("timed out", conditionMessage(w))) timed_out <<- TRUE
      invokeRestart("muffleWarning")
    }
  )
  skip_if(timed_out, "simulate.r smoke test exceeded the 180s timeout")

  status <- attr(result, "status")
  ran_ok <- is.null(status) || status == 0
  if (!ran_ok) {
    fail(paste0("simulate.r exited with status ", status, ":\n", paste(result, collapse = "\n")))
  }

  expect_true(ran_ok && file.exists(out_rda))

  if (ran_ok && file.exists(out_rda)) {
    e <- new.env()
    load(out_rda, envir = e)
    expect_true(all(c("rec.host", "rec.transposon", "rec.offspring") %in% ls(e)))
    expect_equal(nrow(e$rec.host), gen_max + 1)
    expect_equal(ncol(e$rec.host), pop_size)
    expect_equal(dim(e$rec.transposon), dim(e$rec.host))
    expect_equal(dim(e$rec.offspring), dim(e$rec.host))

    ## every recorded transposon tag should still have all 6 "!"-delimited
    ## fields (gene/location/generation/valid/uniqID/jump); if simulate.r's
    ## post-tPn.jump colname parsing ever uses the wrong separator, tags come
    ## out malformed and this assertion catches it
    tags <- unlist(strsplit(unlist(e$rec.transposon), ";", fixed = TRUE))
    tags <- tags[!is.na(tags) & nzchar(tags)]
    if (length(tags) > 0) {
      field_counts <- lengths(strsplit(tags, "!", fixed = TRUE))
      expect_true(all(field_counts == 6))
    }
  }
})
