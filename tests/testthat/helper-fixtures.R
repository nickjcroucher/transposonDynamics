# Shared fixtures for the func.r / simulate.r test suite.
#
# func.r and simulate.r both use bare, hard-coded relative paths
# (e.g. "../raw/tpn-template.csv", "../data/...", source("func.r")) that
# only resolve correctly when R's working directory is the folder the two
# scripts live in, with sibling raw/ and data/ folders one level up:
#
#   <somewhere>/raw/    input.csv, scenario.csv, seed.csv, tpn-template.csv
#   <somewhere>/data/   <accession>.gff (not supplied), simulate.r's *.rda output
#   <somewhere>/<code>/ func.r, simulate.r      <- working directory when run
#
# testthat itself runs with the working directory set to tests/testthat/,
# which does not match that layout. Rather than guess the real repo's
# folder name for <code>, these helpers locate func.r/simulate.r, then copy
# them (and the raw/ CSVs, and any real data/*.gff already in the repo)
# into a throwaway temp directory that mirrors the layout the scripts
# themselves assume. The original files are only ever read, never written
# to or modified.

locate_scripts <- function() {
  candidates <- c(
    ".", "..", "../..", "../../..",
    "../../R", "../R", "../../code", "../code",
    "../../src", "../src", "../../scripts", "../scripts",
    "../../script", "../script"
  )
  for (d in candidates) {
    if (file.exists(file.path(d, "func.r")) && file.exists(file.path(d, "simulate.r"))) {
      return(normalizePath(d))
    }
  }
  # Fall back to a recursive search a couple of levels up from tests/testthat/.
  search_root <- normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE)
  hits <- list.files(search_root, pattern = "^func\\.r$", recursive = TRUE,
                      full.names = TRUE, ignore.case = TRUE)
  if (length(hits) > 0) return(normalizePath(dirname(hits[1])))
  stop(
    "Could not find func.r/simulate.r near the test directory. If this ",
    "repo keeps them somewhere unusual, add that path to the 'candidates' ",
    "vector at the top of tests/testthat/helper-fixtures.R.",
    call. = FALSE
  )
}

# Builds an isolated copy of the expected raw/ + data/ + code/ layout in a
# temp directory and returns its path. Safe to call repeatedly; each call
# gets its own throwaway directory.
build_sandbox <- function() {
  scripts_dir <- locate_scripts()
  raw_dir_src <- file.path(scripts_dir, "..", "raw")
  data_dir_src <- file.path(scripts_dir, "..", "data")

  sandbox <- file.path(
    tempdir(),
    paste0("tpnsim-", format(Sys.time(), "%Y%m%d%H%M%OS3"), "-", sample.int(1e6, 1))
  )
  dir.create(file.path(sandbox, "code"), recursive = TRUE)
  dir.create(file.path(sandbox, "raw"), recursive = TRUE)
  dir.create(file.path(sandbox, "data"), recursive = TRUE)

  file.copy(file.path(scripts_dir, "func.r"), file.path(sandbox, "code", "func.r"))
  file.copy(file.path(scripts_dir, "simulate.r"), file.path(sandbox, "code", "simulate.r"))

  for (f in c("input.csv", "scenario.csv", "seed.csv", "tpn-template.csv")) {
    src <- file.path(raw_dir_src, f)
    if (file.exists(src)) file.copy(src, file.path(sandbox, "raw", f))
  }

  # Pass through a real genome GFF if this repo already has one, so the true
  # end-to-end path (inParams() -> simulate.r) can run for real. Tests that
  # need it call skip_if()/skip_if_not() around its presence rather than
  # failing when it's absent -- see test-func.R / test-simulate.R.
  if (dir.exists(data_dir_src)) {
    gffs <- list.files(data_dir_src, pattern = "\\.gff$", full.names = TRUE)
    for (g in gffs) file.copy(g, file.path(sandbox, "data", basename(g)))
  }

  normalizePath(sandbox)
}

# Sources func.r into its own environment with the working directory set to
# the sandbox's code/ folder, so its bare "../raw/tpn-template.csv" (func.r
# line 11) resolves. Returns the environment holding all of func.r's
# functions and top-level objects (e.g. fenv$tPn.0, fenv$gVar).
source_func_r <- function(sandbox) {
  env <- new.env(parent = globalenv())
  old <- setwd(file.path(sandbox, "code"))
  on.exit(setwd(old), add = TRUE)
  source("func.r", local = env)
  env
}

# Runs simulate.r as a real subprocess (it reads commandArgs() itself, so it
# can't just be source()'d in-process) with the working directory set the
# same way. Returns the combined stdout+stderr lines and the exit status.
run_rscript <- function(sandbox, args = character(0), timeout = 60) {
  old <- setwd(file.path(sandbox, "code"))
  on.exit(setwd(old), add = TRUE)
  out <- suppressWarnings(system2("Rscript", c("simulate.r", args),
                                   stdout = TRUE, stderr = TRUE, timeout = timeout))
  status <- attr(out, "status")
  list(output = out, status = if (is.null(status)) 0L else status)
}

# A small, hand-built gene table matching the structure inParams() would
# otherwise build from a genome GFF (locus_tag/start/end/length/interLength/
# essential/recombination/advantage), so functions downstream of inParams()
# can be tested without the external genome annotation file.
make_gene_df <- function() {
  data.frame(
    locus_tag = c("SPD_0001", "SPD_0002", "SPD_0003"),
    start = c(1, 501, 1101),
    end = c(400, 900, 1400),
    product = c("geneA", "geneB", "geneC"),
    length = c(400, 400, 300),
    interLength = c(100, 200, 150),
    essential = c(TRUE, FALSE, FALSE),
    recombination = c(FALSE, TRUE, FALSE),
    advantage = c(FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
}

# Writes a tiny, fast-running variant of input.csv (small population/
# generation counts) using whichever "ref genome" value is passed in, so a
# full simulate.r run can complete in seconds when a real GFF is available.
make_small_input_csv <- function(path, ref_genome_value) {
  df <- data.frame(
    Type = c(
      "ref genome", "essential genes", "genes for recombination mechanism",
      "genes for fitness advantage", "percentage of fitness benefit with transposon",
      "host organism constant generation number", "host organism constant population size",
      "host genome variation", "host genetic variation", "transposon size in bp",
      "transposon population size per genome mean", "transposon population size per genome sd",
      "transposon population size per genome distribution",
      "percentage transposon perturbation genotoxic",
      "percentage chance transposon perturbation boost",
      "percentage amplitude transposon perturbation boost"
    ),
    Value = c(
      ref_genome_value, "", "", "", "10", "2", "4", "3", "10", "1000",
      "2", "0", "poisson", "0.1", "80", "150"
    ),
    Note = "",
    stringsAsFactors = FALSE
  )
  write.csv(df, path, row.names = FALSE)
  invisible(path)
}
