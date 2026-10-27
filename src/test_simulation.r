#!/bin/env Rscript
# author: ph-u (test suite drafted with Claude)
# script: test_simulation.r
# desc: unit, input-data and end-to-end tests for func.r, func_analysis.r, simulate.r and sim_summary.r
# in: Rscript test_simulation.r                      (run from the folder that holds func.r)
#     or, inside R: testthat::test_file("test_simulation.r")
# out: test report on screen; Rscript exits with status 1 if any test fails or errors
# arg: none; optional environment switch:
#      TPN_SKIP_INTEGRATION=true   skip the end-to-end simulate.r runs (the slow part: 18 short runs)
# date: 20261027
#
# Needs internet access: func.r and func_analysis.r source helper scripts (eta.r, colour.r) from GitHub.
# Input file: simulate.r's default (../data/input.csv) when it exists, otherwise the committed ../raw/input.csv
# (data/ is git-ignored, so CI uses the raw/ copy). Templates always come from ../raw, as func.r reads them.
# The suite only READS the project folders. Everything that writes (mock genomes, simulate.r and sim_summary.r
# outputs) happens in throwaway copies of the project under tempdir(), which R deletes when the session ends.

##### 0. `Rscript test_simulation.r` re-runs this file under testthat #####
if (!identical(Sys.getenv("TESTTHAT"), "true")) {
  if (!requireNamespace("testthat", quietly = TRUE) || utils::packageVersion("testthat") < "3.0.0") {
    stop("testthat >= 3.0.0 is required: install.packages(\"testthat\")", call. = FALSE)
  }
  if (interactive()) {
    stop("Run this file with testthat::test_file(\"test_simulation.r\"), not source().", call. = FALSE)
  }
  self <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  self <- if (length(self) == 1) gsub("~+~", " ", self, fixed = TRUE) else "test_simulation.r"
  testthat::test_file(normalizePath(self), reporter = "progress", stop_on_failure = TRUE)
  quit(save = "no", status = 0)
}

##### 1. Set-up #####
library(testthat)
local_edition(3)
if (!requireNamespace("withr", quietly = TRUE)) stop("Please install the 'withr' package.", call. = FALSE)

code_dir <- normalizePath(getwd())
scripts <- c("func.r", "func_analysis.r", "simulate.r", "sim_summary.r")
if (!all(file.exists(file.path(code_dir, scripts)))) {
  stop("test_simulation.r must sit next to ", paste(scripts, collapse = ", "), " (folder: ", code_dir, ")",
       call. = FALSE)
}
SKIP_INTEGRATION <- tolower(Sys.getenv("TPN_SKIP_INTEGRATION")) %in% c("true", "1", "yes")
skip_integration <- function() skip_if(SKIP_INTEGRATION, "TPN_SKIP_INTEGRATION is set")

# simulate.r's built-in defaults: input file, seed file, seed index, scenario file, scenario index
sim_defaults <- local({
  src <- readLines(file.path(code_dir, "simulate.r"), warn = FALSE)
  ln <- grep("argv\\s*=\\s*c\\(", src, value = TRUE)
  if (length(ln) == 0) NULL else gsub("\"", "", regmatches(ln[1], gregexpr("\"[^\"]*\"", ln[1]))[[1]])
})
if (length(sim_defaults) < 5) stop("could not read the default arguments on simulate.r's argv line", call. = FALSE)
from_code <- function(p) normalizePath(file.path(code_dir, p), mustWork = FALSE)
f_input <- from_code(sim_defaults[1])                                 # simulate.r's default (../data/input.csv)
if (!file.exists(f_input)) f_input <- from_code("../raw/input.csv")   # data/ is git-ignored: use the committed copy
f_seed <- from_code(sim_defaults[2])
f_scen <- from_code(sim_defaults[4])
raw_dir <- from_code("../raw")                                        # func.r always reads the templates from here
data_dir <- from_code("../data")

raw_files <- if (dir.exists(raw_dir)) list.files(raw_dir, full.names = TRUE) else character(0)
f_tpn <- raw_files[grepl("-tpn\\.csv$", raw_files)]
f_host <- raw_files[grepl("-host\\.csv$", raw_files)]
if (!file.exists(f_input) || length(f_tpn) != 1 || length(f_host) != 1) {
  stop("the tests need an input file (", sim_defaults[1], " or ../raw/input.csv) and exactly one *-tpn.csv and ",
       "one *-host.csv in ../raw (found ", length(f_tpn), " and ", length(f_host), ")", call. = FALSE)
}
tpn_tmpl <- utils::read.csv(f_tpn, header = TRUE, na.strings = "")     # read as func.r reads it
host_tmpl <- utils::read.csv(f_host, header = TRUE)

cache <- new.env()                       # results computed once, on first use
real_scenarios <- function() {
  if (is.null(cache$scen)) {
    cache$scen <- utils::read.csv(f_scen, header = TRUE, colClasses = c("character", "integer"), na.strings = "")
  }
  cache$scen
}

##### 2. Throwaway project copies with mock genomes #####
# Gene coordinates from gene lengths and the gap before each gene (a negative gap is an overlap).
# Genes are long (20 kb) so two independent insertions rarely hit the same base by chance.
mock_genome <- function(n, len = 20000, gap = 4000, lead = 5000, trail = 9000, tiny = integer(0), tiny_len = 30,
                        tiny_gap = 10, overlap_at = integer(0), overlap = 4, adjacent_at = integer(0)) {
  lens <- rep(len, n)
  gaps <- rep(gap, n)
  lens[tiny] <- tiny_len
  gaps[tiny] <- tiny_gap
  gaps[1] <- lead
  gaps[adjacent_at] <- 0
  gaps[overlap_at] <- -overlap
  starts <- ends <- numeric(n)
  pos <- 0
  for (i in seq_len(n)) {
    starts[i] <- pos + gaps[i] + 1
    ends[i] <- starts[i] + lens[i] - 1
    pos <- ends[i]
  }
  list(n = n, tags = sprintf("MOCK_%04d", seq_len(n)), starts = starts, ends = ends, G = pos + trail)
}

# GFF with a chromosome region, an untagged repeat_region, gene + CDS lines, and a second (plasmid) region
write_gff <- function(path, g) {
  s <- as.integer(g$starts)
  e <- as.integer(g$ends)
  gene <- sprintf("mockchr\tmock\tgene\t%d\t%d\t.\t+\t.\tID=gene-%s;Name=%s;gbkey=Gene;locus_tag=%s",
                  s, e, g$tags, g$tags, g$tags)
  cds <- sprintf(paste0("mockchr\tmock\tCDS\t%d\t%d\t.\t+\t0\tID=cds-%s;Parent=gene-%s;gbkey=CDS;",
                        "locus_tag=%s;product=mock protein %d"),
                 s, e, g$tags, g$tags, g$tags, seq_len(g$n))
  region <- sprintf(paste0("mockchr\tmock\tregion\t1\t%d\t.\t+\t.\tID=mockchr:1..%d;Dbxref=taxon:0;",
                           "gbkey=Src;mol_type=genomic DNA"),
                    as.integer(g$G), as.integer(g$G))
  rpt <- sprintf("mockchr\tmock\trepeat_region\t%d\t%d\t.\t+\t.\tID=id-rpt1;rpt_type=direct", s[1], s[1] + 99L)
  plasmid <- "mockplasmid\tmock\tregion\t1\t5000\t.\t+\t.\tID=mockplasmid:1..5000;gbkey=Src;mol_type=genomic DNA"
  writeLines(c("##gff-version 3", sprintf("##sequence-region mockchr 1 %d", as.integer(g$G)), region, rpt,
               as.vector(rbind(gene, cds)), plasmid), path)
}

build_sandbox <- function(name, genome, fitness, params, scenarios) {
  root <- file.path(tempdir(), "tpn_tests", name)
  unlink(root, recursive = TRUE)
  sb <- list(root = root, code = file.path(root, "code"), raw = file.path(root, "raw"),
             data = file.path(root, "data"))
  for (d in c(sb$code, sb$raw, sb$data)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  copied <- c(file.copy(file.path(code_dir, scripts), sb$code, overwrite = TRUE),
              file.copy(c(f_tpn, f_host), sb$raw, overwrite = TRUE))
  if (!all(copied)) stop("could not copy the project into ", root)
  if (file.exists(f_seed)) {
    file.copy(f_seed, file.path(sb$raw, "seed.csv"), overwrite = TRUE)
  } else {
    writeLines(c("12345", "67890"), file.path(sb$raw, "seed.csv"))
  }
  pms <- utils::read.csv(f_input, header = TRUE)
  for (key in names(params)) {
    if (!key %in% pms$Type) {
      extra <- pms[1, , drop = FALSE]
      extra[1, ] <- ""
      extra$Type <- key
      pms <- rbind(pms, extra)
    }
    pms$Value[pms$Type == key] <- as.character(params[[key]])
  }
  utils::write.csv(pms, file.path(sb$raw, "input.csv"), row.names = FALSE)
  utils::write.csv(scenarios, file.path(sb$raw, "scenario.csv"), row.names = FALSE)
  write_gff(file.path(sb$data, "mock_genome.gff"), genome)
  utils::write.csv(data.frame(locus_tag = fitness$tag, fitness = fitness$value),
                   file.path(sb$data, "mock_fitness.csv"), row.names = FALSE)
  sb
}

# Main mock genome: 40 genes with a lead and tail gap, a 4-bp overlap (genes 10/11), back-to-back genes (20/21)
# and five 30-bp genes (25-29).
mock_main <- mock_genome(40, tiny = 25:29, overlap_at = 11, adjacent_at = 21)
REC_MAIN <- mock_main$tags[c(25, 26)]     # tiny, so random transposons rarely hit them
fit_main <- local({
  keep <- !seq_len(mock_main$n) %in% c(27, 28)   # left out on purpose: gIn() must give them fitness 0
  list(tag = mock_main$tags[keep],
       value = rep(c(1, 0.9, 0.8, 0.7, 0.6, 0.5), length.out = mock_main$n)[keep])
})
GMAX <- 5
HPOP <- 30
base_params <- list(
  "ref genome" = "mock_genome",
  "gene fitness" = "../data/mock_fitness.csv",
  "genes for recombination mechanism" = paste(REC_MAIN, collapse = ";"),
  "host organism constant generation number" = GMAX,
  "host organism constant population size" = HPOP,
  "transposon population size per genome mean" = 8,
  "full result" = "yes"
)
sb_main <- build_sandbox("main", mock_main, fit_main, base_params,
                         data.frame(transposon = tpn_tmpl$uniqID[1], host = 1L))
SB_ARGV <- c("../raw/input.csv", "../raw/seed.csv", "1", "../raw/scenario.csv", "1")

fenv <- new.env()                          # func.r, sourced inside the main copy; it reads argv on load
assign("argv", SB_ARGV, envir = fenv)
withr::with_dir(sb_main$code, sys.source("func.r", envir = fenv, keep.source = FALSE))
aenv <- new.env()                          # func_analysis.r; its colour.r lands in the global environment
sys.source(file.path(code_dir, "func_analysis.r"), envir = aenv, keep.source = FALSE)
tpl <- fenv$tPn.0                          # transposon template exactly as func.r reads it

##### 3. Helpers for the unit tests #####
mk_tpn <- function(gene, location, size = 1000, generation = 0, copyDir = "both") {
  r <- tpl[1, ]
  r$gene <- gene
  r$location <- location
  r$size <- size
  r$generation <- generation
  r$copyDir <- copyDir
  r$valid <- TRUE
  fenv$tPn.io(r)
}
# Toy gene table laid out [inter 10][gene 20] [inter 5][gene 20] [inter 5][gene 20] = positions 1-30, 31-55, 56-80
toy_gdf <- data.frame(locus_tag = c("0T_0001", "0T_0002", "0T_0003"), start = c(11, 36, 61), end = c(30, 55, 80),
                      length = c(20, 20, 20), interLength = c(10, 5, 5), fitness = 1, recombination = NA)
# Two-copy toy table on a 100-bp circle; gene 1's gap (30) = 10 bases before it + 20 after gene 3
cl_gdf <- data.frame(locus_tag = c(paste0("0T_000", 1:3), paste0("1T_000", 1:3)), start = rep(c(11, 36, 61), 2),
                     end = rep(c(30, 55, 80), 2), length = 20, interLength = rep(c(30, 5, 5), 2),
                     fitness = 1, recombination = NA)
# Five genes per copy: 1 kb apart, then 48 kb and ~1 Mb away
hot <- data.frame(locus_tag = c(sprintf("0H_%04d", 1:5), sprintf("1H_%04d", 1:5)),
                  start = rep(c(1, 1001, 2001, 50001, 1e6), 2))
letters_of <- function(x) unique(strsplit(paste(x, collapse = ""), "")[[1]])
gin_main <- function(seed = 1) withr::with_dir(sb_main$code, withr::with_seed(seed, fenv$gIn()))

##### 4. Unit tests: func.r #####
test_that("xIn() returns numbers for numeric values and text otherwise", {
  ref <- data.frame(Type = c("a", "b", "c", "d"), Value = c("10", "x;y", "1e-3", " 7 "))
  expect_identical(fenv$xIn("a", ref), 10)
  expect_identical(fenv$xIn("b", ref), "x;y")
  expect_equal(fenv$xIn("c", ref), 0.001)
  expect_identical(fenv$xIn("d", ref), 7)
})

test_that("xIn() reads the input file named in argv[1] and stops on unknown names", {
  expect_identical(fenv$xIn("host organism constant population size"), HPOP)
  expect_identical(fenv$xIn("ref genome"), "mock_genome")
  expect_error(fenv$xIn("no such parameter"))
})

test_that("gffClean() reads the one GFF named after 'ref genome' and splits attributes into columns", {
  withr::local_dir(sb_main$code)
  g <- fenv$gffClean()
  expect_identical(names(g)[1:8], c("seqname", "source", "type", "start", "end", "score", "strand", "phase"))
  expect_true(all(c("ID", "Name", "gbkey", "locus_tag", "Parent", "product", "rpt_type") %in% names(g)))
  expect_equal(nrow(g), 3 + 2 * mock_main$n)          # 2 regions + 1 repeat + gene and CDS lines; '#' lines skipped
  expect_equal(sum(is.na(g$locus_tag)), 3)
  expect_identical(g$locus_tag[g$type == "gene"], mock_main$tags)
  expect_equal(g$start[g$type == "gene"], mock_main$starts)
  expect_identical(g$product[g$type == "CDS"][2], "mock protein 2")
})

test_that("gffClean() matches the GFF on the part of 'ref genome' before its first dot", {
  withr::local_dir(sb_main$code)
  orig <- fenv$pMs
  withr::defer(assign("pMs", orig, envir = fenv))
  p <- orig
  p$Value[p$Type == "ref genome"] <- "mock_genome.7"      # a versioned accession still finds mock_genome.gff
  assign("pMs", p, envir = fenv)
  g <- fenv$gffClean()
  expect_identical(g$locus_tag[g$type == "gene"], mock_main$tags)
})

test_that("gffClean() stops unless exactly one GFF matches 'ref genome'", {
  expect_error(fenv$gffClean(fLst = "../raw/seed.csv"), "expected one GFF")
  expect_error(fenv$gffClean(fLst = c("../data/mock_genome.gff", "../data/old_mock_genome.gff")), "expected one GFF")
})

test_that("gIn() builds the two-copy gene table, fitness and recombination flags", {
  gi <- gin_main()
  g <- mock_main
  n <- g$n
  expect_named(gi, c("gene", "transposon.titre", "numGene", "genome"))
  expect_equal(gi$numGene, n)                           # untagged repeat and both region rows dropped
  expect_equal(gi$genome, c(1, g$G))                    # first region (chromosome), not the plasmid
  expect_false(anyNA(gi$gene$locus_tag))
  expect_identical(gi$gene$locus_tag, c(paste0("0", g$tags), paste0("1", g$tags)))
  fit <- fit_main$value[match(g$tags, fit_main$tag)]
  fit[is.na(fit)] <- 0                                  # genes missing from the file count as essential
  expect_equal(gi$gene$fitness, c(fit, fit))
  rec <- g$tags %in% REC_MAIN
  expect_equal(gi$gene$recombination, c(ifelse(rec, 0, NA), ifelse(rec, 1, NA)))
  expect_length(gi$transposon.titre, HPOP)
  expect_true(all(gi$transposon.titre >= 1 & gi$transposon.titre == round(gi$transposon.titre)))
  expect_identical(gi, gin_main())                     # same seed, same result
})

test_that("gIn() lengths and gaps cover the whole circular chromosome", {
  gi <- gin_main()
  g <- mock_main
  n <- g$n
  c0 <- gi$gene[seq_len(n), ]
  len <- g$ends - g$starts + 1
  gap <- c((g$G - g$ends[n]) + (g$starts[1] - 1), g$starts[-1] - g$ends[-n] - 1)
  ov <- which(gap < 0)                                  # overlap: trimmed off the earlier gene, gap 0
  len[ov - 1] <- len[ov - 1] + gap[ov]
  gap[ov] <- 0
  expect_equal(c0$length, len)
  expect_equal(c0$interLength, gap)
  expect_equal(c0$interLength[c(11, 21)], c(0, 0))      # overlapping and back-to-back genes
  expect_equal(c0$length[10], g$ends[10] - g$starts[10] + 1 - 4)
  expect_equal(sum(c0$length + c0$interLength), g$G)
  expect_true(all(gi$gene$interLength >= 0))
  expect_equal(gi$gene$length[n + seq_len(n)], c0$length)
})

test_that("ini.host() returns host.var genomes of nGenes alleles from the allowed alphabet", {
  h <- withr::with_seed(2, fenv$ini.host(nGenes = 12, host.var = 5, gene.var = 3))
  expect_type(h, "character")
  expect_length(h, 5)
  expect_true(all(nchar(h) == 12))
  expect_true(all(letters_of(h) %in% c("A", "B", "C", "a", "b", "c")))
  expect_identical(h, withr::with_seed(2, fenv$ini.host(12, 5, 3)))
  expect_length(withr::with_seed(2, fenv$ini.host(3, 1, 2)), 1)
})

test_that("ini.host() clamps variation to 1-26 per gene (2 to 52 alleles)", {
  expect_true(all(letters_of(withr::with_seed(3, fenv$ini.host(50, 4, 0))) %in% c("A", "a")))
  many <- letters_of(withr::with_seed(3, fenv$ini.host(500, 4, 100)))
  expect_true(all(many %in% c(LETTERS, letters)))
  expect_gt(length(many), 26)
})

test_that("ini.host() takes per-gene variation written as 'a;b;...' and recycles it over the genes", {
  h <- withr::with_seed(4, fenv$ini.host(nGenes = 4, host.var = 200, gene.var = "1;2"))
  expect_true(all(substr(h, 1, 1) %in% c("A", "a")))
  expect_true(all(substr(h, 2, 2) %in% c("A", "B", "a", "b")))
  expect_true(all(substr(h, 3, 3) %in% c("A", "a")))
  expect_true(all(substr(h, 4, 4) %in% c("A", "B", "a", "b")))
  expect_true(any(substr(h, 2, 2) %in% c("B", "b")))
})

test_that("tPn.io() round-trips template rows between data frame and text", {
  rows <- tpl[1:5, ]
  s <- fenv$tPn.io(rows)
  expect_type(s, "character")
  expect_length(s, 1)
  expect_length(strsplit(s, ";", fixed = TRUE)[[1]], 5)
  expect_false(grepl(" ", s, fixed = TRUE))
  back <- fenv$tPn.io(s)
  expect_identical(names(back), names(tpl))
  expect_equal(back, rows, ignore_attr = TRUE)
})

test_that("tPn.io() joins vectors with ';' and drops empty tokens", {
  parts <- strsplit(fenv$tPn.io(tpl[1:3, ]), ";", fixed = TRUE)[[1]]
  expect_equal(nrow(fenv$tPn.io(parts)), 3)
  expect_equal(nrow(fenv$tPn.io(c(parts[1], "", parts[2]))), 2)
})

test_that("tPn.io() keeps rates to about 7 significant digits", {
  r <- tpl[1, ]
  r$jumpRate <- 1 / 3
  r$copyRate <- 0.123456789
  back <- fenv$tPn.io(fenv$tPn.io(r))
  expect_equal(back$jumpRate, 1 / 3, tolerance = 1e-6)
  expect_equal(back$copyRate, 0.123456789, tolerance = 1e-6)
})

test_that("reGeneDF() adds transposon sizes to the gene or intergenic length they sit in", {
  expect_identical(fenv$reGeneDF("", toy_gdf), toy_gdf)
  s <- paste(mk_tpn("g0T_0002", 3, size = 1000), mk_tpn("g0T_0002", 9, size = 500),
             mk_tpn("i0T_0003", 2, size = 100), sep = ";")
  want <- toy_gdf
  want$length[2] <- 20 + 1500
  want$interLength[3] <- 5 + 100
  expect_equal(fenv$reGeneDF(s, toy_gdf), want)
})

test_that("tPn.reloc() returns one transposon string and a gene table grown by its size", {
  t1 <- mk_tpn(NA, 0)
  for (u in c(0.05, 0.3, 0.5, 0.74, 0.9, 0.99999)) {
    r <- fenv$tPn.reloc(num = u, tpn = t1, gEn = 7, gDF = toy_gdf)
    expect_named(r, c("tpn", "gDF"))
    d <- fenv$tPn.io(r$tpn)
    expect_equal(nrow(d), 1)
    expect_equal(d$generation, 7)
    expect_match(d$gene, "^[gi]0T_000[123]$")
    expect_true(d$location >= 1)
    expect_identical(names(r$gDF), names(toy_gdf))
    expect_equal(sum(r$gDF$length + r$gDF$interLength), 80 + 1000)
  }
})

test_that("tPn.reloc() places a position in the region that contains it, including row boundaries", {
  t1 <- mk_tpn(NA, 0)
  at <- function(p) {              # p = 1-based position on the toy table; ceiling(num * 80) == p
    r <- fenv$tPn.reloc(num = (p - 0.5) / 80, tpn = t1, gDF = toy_gdf)
    d <- fenv$tPn.io(r$tpn)
    list(where = c(d$gene, d$location), gDF = r$gDF)
  }
  expect_identical(at(40)$where, c("g0T_0002", "5"))
  expect_identical(at(33)$where, c("i0T_0002", "3"))
  expect_identical(at(4)$where, c("i0T_0001", "4"))
  expect_identical(at(70)$where, c("g0T_0003", "10"))
  expect_identical(at(30)$where, c("g0T_0001", "20"))   # last base of row 1
  expect_identical(at(31)$where, c("i0T_0002", "1"))    # first base of row 2
  expect_identical(at(80)$where, c("g0T_0003", "20"))   # the very last position
  expect_equal(at(40)$gDF$length[2], 20 + 1000)
  expect_equal(at(33)$gDF$interLength[2], 5 + 1000)
})

test_that("tPn.newloc() relocates rows in order and carries each host's updated gene table", {
  seen <- numeric(0)
  stub <- function(num, tpn, gEn = 0, gDF) {
    seen <<- c(seen, sum(gDF$length + gDF$interLength))
    gDF$length[1] <- gDF$length[1] + 1000
    list(tpn = paste0(tpn, "@", gEn), gDF = gDF)
  }
  orig <- fenv$tPn.reloc
  assign("tPn.reloc", stub, envir = fenv)
  withr::defer(assign("tPn.reloc", orig, envir = fenv))
  df <- data.frame(tpn = c("a", "b", "c", "d"), loc = c(0.1, 0.2, 0.3, 0.4), host = c(1, 1, 2, 2), newLoc = NA)
  out <- fenv$tPn.newloc(df, gen = 3, gDF = toy_gdf)
  expect_identical(out$newLoc, c("a@3", "b@3", "c@3", "d@3"))
  expect_equal(seen, c(80, 1080, 80, 1080))           # second row of a host sees the first insertion
})

test_that("tPn.x() invalidates an earlier transposon disrupted by a later insertion in the same region", {
  a <- mk_tpn("g0T_0002", 10, size = 100)               # occupies 10..109
  valid <- function(...) fenv$tPn.io(fenv$tPn.x(paste(..., sep = ";")))$valid
  expect_identical(fenv$tPn.x(""), "")
  expect_identical(valid(a), TRUE)
  expect_identical(valid(a, mk_tpn("g0T_0002", 109, size = 50)), c(FALSE, TRUE))
  expect_identical(valid(a, mk_tpn("g0T_0002", 110, size = 50)), c(TRUE, TRUE))
  expect_identical(valid(a, mk_tpn("g0T_0003", 50, size = 50)), c(TRUE, TRUE))
  expect_identical(valid(sub("!TRUE!", "!FALSE!", a, fixed = TRUE)), TRUE)   # flags are recomputed
})

test_that("tPn.flat() moves notation between chromosome copies without touching other fields", {
  s0 <- mk_tpn("g0T_0002", 10)
  s1 <- fenv$tPn.flat(s0, ori = 0)
  i0 <- mk_tpn("i0T_0003", 4)
  expect_identical(s1, sub("^g0", "g1", s0))
  expect_identical(fenv$tPn.flat(s1, ori = 1), s0)
  expect_identical(fenv$tPn.flat(s0, ori = 2), s1)     # ori is taken modulo 2
  expect_identical(fenv$tPn.flat(c(s0, i0), 0), c(s1, sub("^i0", "i1", i0)))
  expect_identical(fenv$tPn.flat(s1, ori = 0), s1)     # nothing on copy 0 to move
})

test_that("h1.mod() leaves 'fixed' rates unchanged and always returns numbers", {
  expect_identical(fenv$h1.mod(0.1, "fixed"), 0.1)
  expect_identical(fenv$h1.mod(c("0.1", " 1e-04"), "fixed"), c(0.1, 1e-04))   # text, as apply() passes it
})

test_that("h1.mod() divides by copy number under 'charlesworth'", {
  expect_equal(fenv$h1.mod(0.1, "charlesworth", numTpn = 4), 0.025)
  expect_equal(fenv$h1.mod(c(0.1, 0.2), "charlesworth", numTpn = c(2, 4)), c(0.05, 0.05))
  expect_equal(fenv$h1.mod(" 0.1", "charlesworth", numTpn = " 5"), 0.02)
  expect_error(fenv$h1.mod(c(0.1, 0.2, 0.3), "charlesworth", numTpn = c(1, 2)), "Length of numTpn")
})

test_that("h1.mod() 'evolving' multiplies each rate by its own |N(1, 0.01)| draw", {
  r <- c(0.1, 0.01, 0)
  got <- withr::with_seed(11, fenv$h1.mod(r, "evolving"))
  want <- withr::with_seed(11, r * abs(rnorm(3, mean = 1, sd = 0.01)))
  expect_identical(got, want)
  expect_identical(got[3], 0)
  expect_true(all(got >= 0))
  steps <- withr::with_seed(12, fenv$h1.mod(rep(1, 5000), "evolving"))
  expect_true(abs(stats::sd(steps) - 0.01) < 0.001, info = "a 1% step whatever the rate")
})

test_that("tPn.comLoc() calls a copy 'terminus' when it lands farther from the origin", {
  cmp <- function(old, new) {
    fenv$tPn.comLoc(c(mk_tpn(old[1], as.numeric(old[2]), copyDir = "terminus"),
                      mk_tpn(new[1], as.numeric(new[2]), copyDir = "terminus")), G = 100, gDF = cl_gdf)
  }
  # positions on the 100-bp circle: gene start - 1 + location, minus the gap length for intergenic labels
  expect_identical(cmp(c("g0T_0002", 5), c("g0T_0002", 10)), "terminus")  # 40 -> 45
  expect_identical(cmp(c("g0T_0002", 5), c("g0T_0002", 1)), "origin")     # 40 -> 36
  expect_identical(cmp(c("g0T_0002", 5), c("g0T_0003", 10)), "origin")    # 40 -> 70: past the terminus, 30 from origin
  expect_identical(cmp(c("g0T_0003", 10), c("g0T_0003", 1)), "terminus")  # 70 -> 61: 30 -> 39 from the origin
  expect_identical(cmp(c("g0T_0002", 5), c("i0T_0002", 1)), "origin")     # 40 -> 31
  expect_identical(cmp(c("g0T_0002", 5), c("i0T_0001", 5)), "origin")     # 40 -> -15, wraps to 85: 15 from origin
  expect_identical(cmp(c("g0T_0002", 5), c("g1T_0002", 10)), "terminus")  # copy digit does not change the position
  expect_identical(cmp(c("g0T_0002", 5), c("g0T_0002", 5)), "terminus")   # same spot: returns copyDir
  expect_error(fenv$tPn.comLoc(mk_tpn("g0T_0002", 5), G = 100, gDF = cl_gdf), "Length-two")
})

test_that("gpSel() grows numbered blocks from the seeds (x == 0) through neighbours with x <= lim", {
  ex <- c(0, .5, .6, .2, .7, .4, .3, .5, .2, .5, .3, 0, .5, .6, .2, .7, .4, .3, .5, .2, .5)   # example in func.r
  expect_equal(fenv$gpSel(ex, rep(.5, 21)), c(2, 2, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 2, 2, 2, 2, 2))
  expect_equal(fenv$gpSel(c(0, .9, 0, .9), rep(.5, 4)), c(1, 0, 2, 0))           # two separate blocks
  expect_equal(fenv$gpSel(c(0, .9, .1, .1, .9), rep(.5, 5)), c(1, 0, 0, 0, 0))   # qualifying genes cut off from seeds
})

test_that("gpSel() treats the chromosome as circular in both directions", {
  expect_equal(fenv$gpSel(c(0, .9, .9, .4), rep(.5, 4)), c(1, 0, 0, 1))   # seed on gene 1 reaches the last gene
  expect_equal(fenv$gpSel(c(.4, .9, .9, 0), rep(.5, 4)), c(1, 0, 0, 1))   # seed on the last gene reaches gene 1
})

test_that("gpSel() edge cases: no seed, every gene qualifies, a single gene", {
  expect_equal(fenv$gpSel(c(.1, .2, .3), rep(.5, 3)), c(0, 0, 0))
  expect_equal(fenv$gpSel(c(0, .1, .1, .1), rep(.5, 4)), c(1, 1, 1, 1))
  expect_equal(fenv$gpSel(0, .5), 1)
})

test_that("gRc.prob() returns each seed gene and its contiguous neighbours, numbered by block", {
  for (s in 1:25) {
    out <- withr::with_seed(s, fenv$gRc.prob("0H_0002", g.hotspot = hot))
    expect_named(out, c("locus_tag", "contiguous"))
    expect_true("0H_0002" %in% out$locus_tag)
    expect_true(all(startsWith(out$locus_tag, "0")))              # stays on the seed's copy
    expect_false(any(c("0H_0004", "0H_0005") %in% out$locus_tag)) # 48 kb or more from the seed
    expect_true(all(out$contiguous == 1))                         # one seed, one block
  }
  expect_equal(nrow(withr::with_seed(1, fenv$gRc.prob(character(0), g.hotspot = hot))), 0)
})

test_that("gRc.prob() joins two seeds through a gene 1 kb from both", {
  two <- withr::with_seed(2, fenv$gRc.prob(c("0H_0001", "0H_0003"), g.hotspot = hot))
  expect_setequal(two$locus_tag, c("0H_0001", "0H_0002", "0H_0003"))   # 0H_0002: threshold 0.5 + 0.5 = 1
  expect_true(all(two$contiguous == 1))
})

test_that("gRc.prob() numbers copy-1 blocks after the copy-0 blocks", {
  out <- withr::with_seed(3, fenv$gRc.prob(c("0H_0002", "1H_0004"), g.hotspot = hot))
  b0 <- out$contiguous[startsWith(out$locus_tag, "0")]
  b1 <- out$contiguous[startsWith(out$locus_tag, "1")]
  expect_true("1H_0004" %in% out$locus_tag)
  expect_true(length(b0) > 0 && length(b1) > 0)
  expect_true(min(b1) > max(b0))
})

test_that("gRc.prob() takes a neighbour 1 kb away about half of the time", {
  hits <- vapply(1:2000, function(s) {
    withr::with_seed(s, "0H_0001" %in% fenv$gRc.prob("0H_0002", g.hotspot = hot)$locus_tag)
  }, logical(1))
  expect_true(abs(mean(hits) - 0.5) < 0.05, info = paste("observed", mean(hits)))
})

test_that("gRc.prob() blocks on the mock genome are unbroken runs on one chromosome copy", {
  gi <- gin_main()
  g <- gi$gene
  n <- gi$numGene
  for (s in 1:10) {
    seeds <- withr::with_seed(s, sample(g$locus_tag, 6))
    out <- withr::with_seed(s, fenv$gRc.prob(seeds, g.hotspot = g))
    expect_true(all(seeds %in% out$locus_tag))
    for (b in unique(out$contiguous)) {
      tags <- out$locus_tag[out$contiguous == b]
      cp <- unique(substr(tags, 1, 1))
      expect_length(cp, 1)                                         # one copy per block
      idx <- sort(match(tags, g$locus_tag[substr(g$locus_tag, 1, 1) == cp[1]]))
      steps <- diff(c(idx, idx[1] + n))                           # gaps around the circle
      expect_equal(sum(steps > 1), as.integer(length(idx) < n), info = paste(tags, collapse = " "))
    }
  }
})

##### 5. Unit tests: func_analysis.r #####
test_that("the colour palette sim_summary.r plots with (cBp) is available after sourcing func_analysis.r", {
  expect_true(exists("cBp", envir = aenv))                  # colour.r is sourced into the global environment
  pal <- get("cBp", envir = aenv)
  expect_gte(length(pal), 1)
  expect_error(grDevices::col2rgb(pal, alpha = TRUE), NA)
})

test_that("fAmily() follows ancestry through copy-0 (1..n) and copy-1 (n+1..2n) parent indices", {
  tr <- rbind(c("1;1", "2;2", "3;3"),
              c("4;4", "4;4", "2;6"))
  out <- aenv$fAmily(tr)
  expect_equal(unlist(out[1, ], use.names = FALSE), c(1, 1, 1))
  expect_equal(unlist(out[2, ], use.names = FALSE), c(2, 0.5, 0.5))
  n <- 40                                                # a population size other than 30
  big <- rbind(paste0(seq_len(n), ";", seq_len(n)), rep("45;45", n))
  expect_equal(unlist(aenv$fAmily(big)[2, ], use.names = FALSE), replace(numeric(n), 5, n))
})

test_that("fAmily() handles a single generation", {
  out <- aenv$fAmily(rbind(c("1;2", "2;2")))
  expect_equal(dim(out), c(1, 2))
  expect_equal(unlist(out[1, ], use.names = FALSE), c(0.5, 1.5))
})

test_that("fAmily() conserves total ancestry in every generation", {
  withr::local_seed(9)
  n <- 12
  g <- 6
  pick <- function(k, m) paste(sample(k, m, TRUE), sample(k, m, TRUE), sep = ";")
  tr <- rbind(pick(3, n), matrix(pick(2 * n, n * (g - 1)), nrow = g - 1))
  expect_equal(unname(rowSums(aenv$fAmily(tr))), rep(n, g))
})

test_that("fAmily(state = 'diploid') counts each genome's founding template, including lost lineages", {
  tr <- rbind(c("1;2", "3;1", "2;2"),       # genomes: copy 0 = 1, 3, 2; copy 1 = 2, 1, 2
              c("1;4", "5;1", "4;4"))       # parents carry only templates 1 and 2: template 3 is lost
  out <- aenv$fAmily(tr, "diploid")
  expect_equal(dim(out), c(2, 3))
  expect_equal(unlist(out[1, ], use.names = FALSE), c(1, 1.5, 0.5))   # genome counts / 2
  expect_equal(unlist(out[2, ], use.names = FALSE), c(1.5, 1.5, 0))
  one <- aenv$fAmily(tr[1, , drop = FALSE], "diploid")                # a single generation
  expect_equal(dim(one), c(1, 3))
})

test_that("fAmily(state = 'diploid') conserves ancestry: each generation sums to the population size", {
  withr::local_seed(19)
  n <- 12
  g <- 6
  pick <- function(k, m) paste(sample(k, m, TRUE), sample(k, m, TRUE), sep = ";")
  first <- c("1;2", "3;1", "2;3", pick(3, n - 3))                    # every founding template is used
  tr <- rbind(first, matrix(pick(2 * n, n * (g - 1)), nrow = g - 1))
  expect_equal(unname(rowSums(aenv$fAmily(tr, "diploid"))), rep(n, g))
})

test_that("fAmily(state = 'diploid') gives one column per founding template actually used", {
  out <- aenv$fAmily(rbind(c("2;3", "3;2", "3;3")), "diploid")       # template 1 never drawn
  expect_equal(unlist(out[1, ], use.names = FALSE), c(1, 2))         # columns = templates 2 and 3
})

##### 6. Input data used by simulate.r's defaults #####
test_that("the input file defines every parameter the scripts read, exactly once", {
  pms <- utils::read.csv(f_input, header = TRUE)
  expect_true(all(c("Type", "Value") %in% names(pms)))
  src <- unlist(lapply(file.path(code_dir, scripts), readLines, warn = FALSE))
  hits <- unlist(regmatches(src, gregexpr('xIn\\("[^"]+"', src)))
  keys <- unique(sub('^xIn\\("', "", sub('"$', "", hits)))
  expect_gt(length(keys), 0)
  n <- vapply(keys, function(k) sum(pms$Type == k), integer(1))
  expect_true(all(n == 1), info = paste("missing or duplicated:", paste(keys[n != 1], collapse = " | ")))
})

test_that("input values are in the ranges the code assumes", {
  pms <- utils::read.csv(f_input, header = TRUE)
  val <- function(k) as.character(pms$Value[pms$Type == k][1])
  num <- function(k) suppressWarnings(as.numeric(val(k)))
  whole <- function(x) length(x) > 0 && all(is.finite(x) & x == round(x))
  gens <- num("host organism constant generation number")
  pop <- num("host organism constant population size")
  gv <- suppressWarnings(as.numeric(strsplit(val("host genetic variation"), ";", fixed = TRUE)[[1]]))
  expect_true(whole(gens) && gens >= 1, info = "generation number")
  expect_true(whole(pop) && pop >= 3, info = "population size (recombination samples a donor among the other hosts)")
  expect_true(whole(num("host genome variation")) && num("host genome variation") >= 1)
  expect_true(whole(gv) && all(gv >= 1 & gv <= 26), info = "host genetic variation: 1-26, optionally per gene with ';'")
  expect_true(isTRUE(num("transposon population size per genome mean") > 0))
  for (k in c("percentage transposon perturbation genotoxic", "percentage chance transposon perturbation boost")) {
    expect_true(isTRUE(num(k) >= 0 && num(k) <= 100), info = k)
  }
  expect_true(isTRUE(num("percentage amplitude transposon perturbation boost") >= 0))
  expect_true(tolower(val("full result")) %in% c("yes", "no", "true", "false"))
})

test_that("the gene fitness file named in the input file exists", {
  skip_if_not(dir.exists(data_dir), "../data not found")
  pms <- utils::read.csv(f_input, header = TRUE)
  fit <- pms$Value[pms$Type == "gene fitness"][1]
  expect_true(file.exists(file.path(code_dir, fit)), info = fit)
})

test_that("transposon template keeps the column order simulate.r reads by position", {
  # simulate.r reads x[11], x[13], x[15]-x[17] (rate rules) and x[1], x[23] (copy direction) by position
  expect_identical(names(tpn_tmpl), c("gene", "location", "generation", "valid", "uniqID", "size",
                                      "jumpRate", "jumpH1", "copyRate", "copyH1", "copyDir"))
})

test_that("transposon template values are valid", {
  t <- tpn_tmpl
  expect_true(all(t$valid %in% TRUE))
  expect_true(all(t$size > 0 & t$size == round(t$size)))
  for (r in c("jumpRate", "copyRate")) expect_true(all(t[[r]] >= 0 & t[[r]] <= 1), info = r)
  for (h in c("jumpH1", "copyH1")) expect_true(all(t[[h]] %in% c("fixed", "charlesworth", "evolving")), info = h)
  expect_true(all(t$copyDir %in% c("both", "terminus", "origin")))
})

test_that("transposon IDs are unique and safe inside transposon strings", {
  ids <- utils::read.csv(f_tpn, header = TRUE, colClasses = "character", na.strings = character(0))$uniqID
  expect_false(anyDuplicated(ids) > 0)
  # tPn.io() re-reads single transposons, so an ID must not look like NA, a logical or a number
  convertible <- ids[vapply(ids, function(u) !is.character(utils::type.convert(u, as.is = TRUE)), logical(1))]
  expect_true(length(convertible) == 0,
              info = paste("tPn.io() would read these as NA/logical/number:", paste(convertible, collapse = ", ")))
  expect_false(any(grepl("[!; ]", ids)), info = "IDs must not contain '!', ';' or spaces")
  expect_false(any(grepl("[gi][01]", ids)), info = "tPn.flat() rewrites g0/g1/i0/i1 anywhere in a transposon string")
})

test_that("host template has the expected columns and values", {
  h <- host_tmpl
  expect_true(all(c("recom", "recomH1", "cell", "homologousAutoRecom", "hAReomH1", "transposonEffect",
                    "genotoxic") %in% names(h)))
  expect_true(all(h$recom >= 0 & h$recom <= 1))
  expect_true(all(h$recomH1 %in% c("switch", "homeostatic")))
  expect_true(all(h$cell %in% c("haploid", "diploid")))
  expect_true(all(h$homologousAutoRecom >= 0 & h$homologousAutoRecom <= 1))
  expect_true(all(h$hAReomH1 %in% c("fixed", "evolving")))
  expect_true(is.logical(h$transposonEffect) && !anyNA(h$transposonEffect))
  expect_true(all(h$genotoxic >= 0 & h$genotoxic == round(h$genotoxic)))
})

test_that("scenario file points at existing template rows", {
  skip_if_not(file.exists(f_scen), "scenario file not found")
  sc <- real_scenarios()
  expect_true(all(c("transposon", "host") %in% names(sc)))
  expect_true(all(!is.na(sc$host) & sc$host >= 1 & sc$host <= nrow(host_tmpl)))
  ids <- unique(unlist(strsplit(unique(sc$transposon), ";", fixed = TRUE)))
  # same rule as simulate.r: tPn.0$uniqID %in% strsplit(sCene$transposon, ";")[[1]]
  n_rows <- vapply(ids, function(id) sum(tpn_tmpl$uniqID %in% id), integer(1))
  expect_true(all(n_rows == 1),
              info = paste("unmatched or ambiguous IDs:", paste(utils::head(ids[n_rows != 1], 10), collapse = ", ")))
})

test_that("simulate.r's default seed and scenario exist, and its skip-read returns that row", {
  skip_if_not(file.exists(f_seed) && file.exists(f_scen), "seed or scenario file not found")
  seeds <- utils::read.csv(f_seed, header = FALSE)[, 1]
  k_seed <- suppressWarnings(as.numeric(sim_defaults[3]))
  k <- suppressWarnings(as.numeric(sim_defaults[5]))
  expect_true(isTRUE(k_seed >= 1 && k_seed <= length(seeds)), info = sim_defaults[3])
  sc <- real_scenarios()
  expect_true(isTRUE(k >= 1 && k <= nrow(sc)), info = sim_defaults[5])
  skip_if_not(isTRUE(k >= 1 && k <= nrow(sc)))
  one <- utils::read.csv(f_scen, header = FALSE, skip = k, nrows = 1, col.names = c("transposon", "host"),
                         colClasses = c("character", "integer"), na.strings = "")
  expect_equal(one, sc[k, ], ignore_attr = TRUE)
})

test_that("seed file holds distinct whole-number seeds", {
  skip_if_not(file.exists(f_seed), "seed file not found")
  s <- utils::read.csv(f_seed, header = FALSE)[, 1]
  expect_gt(length(s), 0)
  expect_true(is.numeric(s) && !anyNA(s) && all(s == round(s)))
  expect_false(anyDuplicated(s) > 0)
})

##### 7. Real genome in ../data (skipped when absent) #####
real_ref <- local({
  pms <- utils::read.csv(f_input, header = TRUE)
  as.character(pms$Value[pms$Type == "ref genome"][1])
})
real_key <- strsplit(real_ref, "[.]")[[1]][1]       # gffClean() matches on the part before the first dot
real_gffs <- if (dir.exists(data_dir)) list.files(data_dir, pattern = "gff$") else character(0)
skip_no_genome <- function() skip_if(length(real_gffs) == 0, "no *.gff in ../data")
real_env <- function() {
  if (is.null(cache$renv)) {
    e <- new.env()
    real_argv <- sim_defaults
    real_argv[1] <- f_input
    assign("argv", real_argv, envir = e)
    withr::with_dir(code_dir, sys.source("func.r", envir = e, keep.source = FALSE))
    cache$renv <- e
  }
  cache$renv
}
real_gff <- function() {                 # parsing a full GFF is slow, so do it once
  if (is.null(cache$rgff)) cache$rgff <- withr::with_dir(code_dir, real_env()$gffClean())
  cache$rgff
}

test_that("real genome: exactly one GFF in ../data is named after 'ref genome'", {
  skip_no_genome()
  expect_equal(sum(grepl(real_key, real_gffs, fixed = TRUE)), 1,
               info = paste("GFF files:", paste(real_gffs, collapse = ", ")))
})

test_that("real genome: a region row, and unique, ordered locus tags", {
  skip_no_genome()
  g <- real_gff()
  expect_true(any(g$type == "region"), info = "gIn() reads the genome size from the first region row")
  u <- unique(g[which(g$type != "region" & !is.na(g$locus_tag)), c("locus_tag", "start", "end")])
  dup <- unique(u$locus_tag[duplicated(u$locus_tag)])
  expect_true(length(dup) == 0,
              info = paste("locus tags with several coordinates:", paste(utils::head(dup), collapse = ", ")))
  expect_false(is.unsorted(u$start), info = "gIn() assumes genes are sorted by start")
})

test_that("real genome: recombination genes listed in the input file exist", {
  skip_no_genome()
  rec <- strsplit(real_env()$xIn("genes for recombination mechanism"), ";", fixed = TRUE)[[1]]
  missing <- setdiff(rec, real_gff()$locus_tag)
  expect_true(length(missing) == 0, info = paste("not in the GFF:", paste(missing, collapse = ", ")))
})

test_that("real genome: the gene fitness file is usable", {
  skip_no_genome()
  f <- file.path(code_dir, real_env()$xIn("gene fitness"))
  skip_if_not(file.exists(f), "gene fitness file not found")
  fit <- utils::read.csv(f, header = TRUE)
  expect_true(all(c("locus_tag", "fitness") %in% names(fit)))
  expect_true(is.numeric(fit$fitness))
  expect_true(all(fit$fitness[!is.na(fit$fitness)] >= 0))
  tags <- unique(real_gff()$locus_tag)
  expect_gt(mean(tags[!is.na(tags)] %in% fit$locus_tag), 0)
})

test_that("real genome: gIn() builds a consistent gene table", {
  skip_no_genome()
  e <- real_env()
  orig <- e$gffClean
  assign("gffClean", function(...) real_gff(), envir = e)     # reuse the parsed GFF
  withr::defer(assign("gffClean", orig, envir = e))
  gi <- withr::with_dir(code_dir, withr::with_seed(1, e$gIn()))
  expect_gt(gi$numGene, 0)
  expect_equal(nrow(gi$gene), 2 * gi$numGene)
  expect_false(anyNA(gi$gene$locus_tag))
  expect_true(all(gi$gene$length > 0))
  expect_true(all(gi$gene$interLength >= 0))
  expect_length(gi$genome, 2)
})

##### 8. End-to-end runs of simulate.r + sim_summary.r (in the throwaway copies) #####
pick_row <- function(tmpl, crit, id = FALSE) {
  ok <- if (id) !is.na(tmpl$uniqID) else rep(TRUE, nrow(tmpl))
  for (nm in names(crit)) ok <- ok & !is.na(tmpl[[nm]]) & tmpl[[nm]] == crit[[nm]]
  hit <- which(ok)[1]
  if (is.na(hit)) NA else if (id) tmpl$uniqID[hit] else hit
}

# One row per scenario: c(transposon uniqID, host template row), picked by property from your templates.
scenario_rows <- function() {
  FX <- "fixed"
  tpn <- function(j, cp, jh = FX, ch = FX, dir = "both") {
    pick_row(tpn_tmpl, list(jumpRate = j, copyRate = cp, jumpH1 = jh, copyH1 = ch, copyDir = dir), id = TRUE)
  }
  host <- function(...) pick_row(host_tmpl, list(...))
  hap <- host(recom = 0.1, recomH1 = "switch", cell = "haploid", hAReomH1 = FX, genotoxic = 0, homologousAutoRecom = 0.9)
  hap_rec <- host(recom = 0.1, recomH1 = "switch", cell = "haploid", hAReomH1 = FX, genotoxic = 0,
                  homologousAutoRecom = 0)
  dip_rec <- host(recom = 0.1, recomH1 = "switch", cell = "diploid", hAReomH1 = FX, genotoxic = 0)
  def <- c(NA, NA)
  k <- suppressWarnings(as.numeric(sim_defaults[5]))
  if (file.exists(f_scen)) {
    sc <- real_scenarios()
    if (isTRUE(k >= 1 && k <= nrow(sc))) def <- c(sc$transposon[k], sc$host[k])
  }
  list(
    default = def,
    fixed = c(tpn(0.1, 0.1), hap),
    diploid = c(tpn(0.1, 0.1), host(recom = 0.1, recomH1 = "switch", cell = "diploid", hAReomH1 = FX,
                                    genotoxic = 0, homologousAutoRecom = 0.5)),
    h1_rules = c(tpn(0.1, 0.1, jh = "charlesworth", ch = "evolving"), hap),
    genotoxic = c(tpn(0.1, 0.1), host(recom = 0.1, recomH1 = "switch", cell = "haploid", hAReomH1 = FX, genotoxic = 4)),
    no_jumps = c(tpn(0, 0.1), hap),
    copy_dir = c(tpn(0.1, 0.1, dir = "terminus"), hap),
    copy_dir_copies_only = c(tpn(0, 0.1, dir = "terminus"), hap),
    copy_dir_jumps_only = c(tpn(0.1, 0, dir = "terminus"), hap),
    recom_tiny = c(tpn(0.1, 0.1), host(recom = 0.0001, recomH1 = "switch", cell = "haploid", hAReomH1 = FX,
                                       genotoxic = 0)),
    recomb_only = c(tpn(0.1, 0.1), hap_rec),
    repair_only = c(tpn(0.1, 0.1), host(recom = 0, cell = "haploid", hAReomH1 = FX, genotoxic = 0,
                                        homologousAutoRecom = 0.9)),
    inert = c(tpn(0, 0), host(recom = 0, cell = "diploid", hAReomH1 = FX, genotoxic = 0)),
    hgt = c(tpn(0, 0), hap_rec),            # inert transposons: only inheritance and recombination change them
    conversion = c(tpn(0, 0), dip_rec)
  )
}

integration <- function() {             # writes the main copy's scenario.csv on first use
  if (is.null(cache$integ)) {
    rows <- scenario_rows()
    ok <- vapply(rows, function(r) !anyNA(r), logical(1))
    if (!any(ok)) stop("no end-to-end scenario could be matched to the templates")
    fill <- rows[[which(ok)[1]]]
    tab <- data.frame(transposon = vapply(rows, function(r) as.character(if (anyNA(r)) fill[1] else r[1]), character(1)),
                      host = vapply(rows, function(r) as.integer(if (anyNA(r)) fill[2] else r[2]), integer(1)),
                      row.names = NULL)
    utils::write.csv(tab, file.path(sb_main$raw, "scenario.csv"), row.names = FALSE)
    cache$integ <- list(index = stats::setNames(seq_along(rows), names(rows)), ok = ok, tab = tab)
  }
  cache$integ
}

sb_extra <- function(name) {            # "sparse" and "lethal" copies, built on first use
  key <- paste0("sb_", name)
  if (is.null(cache[[key]])) {
    FX <- "fixed"
    aa <- pick_row(tpn_tmpl, list(jumpRate = 0.1, copyRate = 0.1, jumpH1 = FX, copyH1 = FX, copyDir = "both"),
                   id = TRUE)
    host <- pick_row(host_tmpl, list(recom = 0, cell = "haploid", hAReomH1 = FX, genotoxic = 0,
                                     homologousAutoRecom = 0))
    if (anyNA(c(aa, host))) skip("no template rows for this scenario")
    p <- base_params
    if (name == "sparse") {
      # mean 0.01: rpois() is redrawn until >= 1, so nearly every host starts with exactly one transposon,
      # and chromosome copies without transposons appear quickly
      p[["transposon population size per genome mean"]] <- 0.01
      p[["host organism constant generation number"]] <- 10
      g <- mock_main
      fit <- fit_main
    } else {
      # lethal: genes back to back (1-bp gaps) and all essential, so generation 1 leaves no viable host
      g <- mock_genome(20, gap = 1, lead = 0, trail = 0, tiny = 1, tiny_len = 3)
      fit <- list(tag = g$tags, value = rep(0, g$n))
      p[["genes for recombination mechanism"]] <- g$tags[10]
    }
    cache[[key]] <- build_sandbox(name, g, fit, p, data.frame(transposon = aa, host = as.integer(host)))
  }
  cache[[key]]
}

run_cache <- new.env()
run_sim <- function(sb, k, seed_index = 1L, input = "../raw/input.csv") {
  key <- paste(basename(sb$root), k, seed_index, input, sep = "_")
  if (is.null(run_cache[[key]])) {
    rscript <- file.path(R.home("bin"), "Rscript")
    args <- c("simulate.r", input, "../raw/seed.csv", as.character(seed_index), "../raw/scenario.csv", as.character(k))
    out <- withr::with_envvar(c(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep)),
                              withr::with_dir(sb$code, suppressWarnings(
                                system2(rscript, args, stdout = TRUE, stderr = TRUE, timeout = 3600))))
    status <- attr(out, "status")
    run_cache[[key]] <- list(
      status = if (is.null(status)) 0L else status,
      log = out,
      tail = paste(c("--- last lines of the simulate.r log ---", utils::tail(out, 30)), collapse = "\n"),
      sim = file.path(sb$data, sprintf("tPn--%s_%s.rda", seed_index, k)),
      pdf = file.path(sb$data, sprintf("ana--%s-%s.pdf", seed_index, k)),
      ana = file.path(sb$data, sprintf("ana--%s_%s.rda", seed_index, k)))
  }
  run_cache[[key]]
}

load_rda <- function(path) {
  e <- new.env()
  load(path, envir = e)
  e
}
# Split recorded transposon strings into a character matrix of the 11 template fields.
tpn_fields <- function(cells) {
  tok <- unlist(strsplit(cells[!is.na(cells)], ";", fixed = TRUE))
  tok <- tok[nzchar(tok)]
  f <- strsplit(tok, "!", fixed = TRUE)
  list(n = lengths(f), m = if (length(f)) do.call(rbind, f[lengths(f) == 11]) else NULL)
}
# Sorted "region+locus@location" keys per chromosome copy (copy digit removed) for one host record.
chrom_keys <- function(cell) {
  out <- list(character(0), character(0))
  m <- if (is.na(cell) || !nzchar(cell)) NULL else tpn_fields(cell)$m
  if (is.null(m)) return(out)
  key <- paste0(substr(m[, 1], 1, 1), substring(m[, 1], 3), "@", as.numeric(m[, 2]))
  copy <- substr(m[, 1], 2, 2)
  list(sort(key[copy == "0"]), sort(key[copy == "1"]))
}
# Identity keys per chromosome copy (region, locus, location, generation, transposon ID; copy digit removed).
copy_keys <- function(cell) {
  m <- if (is.na(cell) || !nzchar(cell)) NULL else tpn_fields(cell)$m
  if (is.null(m)) return(list(character(0), character(0)))
  key <- paste(substr(m[, 1], 1, 1), substring(m[, 1], 3), as.numeric(m[, 2]), m[, 3], m[, 5], sep = "|")
  copy <- substr(m[, 1], 2, 2)
  list(key[copy == "0"], key[copy == "1"])
}
genome_keys <- function(tr_row, g) {    # keys of parent genome g (1..n copy 0, n+1..2n copy 1)
  copy_keys(tr_row[(g - 1) %% length(tr_row) + 1])[[(g - 1) %/% length(tr_row) + 1]]
}
# Host records holding the same transposon twice (identical in every field except 'valid').
same_spot_duplicates <- function(tr) {
  bad <- character(0)
  for (r in seq_len(nrow(tr))) {
    for (j in seq_len(ncol(tr))) {
      cell <- tr[r, j]
      if (is.na(cell) || !nzchar(cell)) next
      tok <- strsplit(cell, ";", fixed = TRUE)[[1]]
      tok <- tok[nzchar(tok)]
      key <- vapply(strsplit(tok, "!", fixed = TRUE), function(f) paste(f[-4], collapse = "!"), character(1))
      if (anyDuplicated(key) > 0) {
        bad <- c(bad, sprintf("generation %d host %d: %s", r - 1, j, key[duplicated(key)][1]))
      }
    }
  }
  bad
}

expect_sim_outputs <- function(run, uniq, cell, recom_rate = NULL, gmax = GMAX, hpop = HPOP) {
  expect_equal(run$status, 0, info = run$tail)
  expect_true(any(grepl("simulation completed", run$log, fixed = TRUE)), info = run$tail)
  if (run$status != 0) return(invisible(NULL))
  for (f in c(run$sim, run$pdf, run$ana)) expect_true(file.exists(f), info = basename(f))
  if (!file.exists(run$sim)) return(invisible(NULL))
  e <- load_rda(run$sim)
  recs <- c("rec.host", "rec.transposon", "rec.offspring", "rec.gRecom", "rec.tpnCount")
  expect_true(all(c("hOst", "sim.df", recs) %in% ls(e)))
  for (m in recs) expect_equal(dim(e[[m]]), c(gmax + 1, hpop), info = m)
  expect_named(e$sim.df, c("familyTree", "hARecomRate", "transposon", "tpnFail", "host"))
  expect_equal(nrow(e$sim.df), hpop)

  # familyTree "a;b": generation 0 points at genome templates, later generations at the 2 x hpop parent genomes
  ft <- e$rec.offspring
  expect_true(all(grepl("^[0-9]+;[0-9]+$", ft)), info = "familyTree entries look like 'a;b'")
  pairs <- matrix(as.integer(unlist(strsplit(as.vector(ft), ";", fixed = TRUE))), nrow = 2)
  first <- rep(seq_len(nrow(ft)), times = ncol(ft)) == 1
  expect_true(all(pairs[, first] >= 1 & pairs[, first] <= length(e$hOst)))
  expect_true(all(pairs[, !first] >= 1 & pairs[, !first] <= 2 * hpop))
  if (identical(cell, "haploid")) {
    expect_true(all(pairs[1, ] == pairs[2, ]), info = "haploid offspring copy one parent genome")
  }

  # transposon records parse into the 11 template fields and stay consistent
  tr <- e$rec.transposon
  tf <- tpn_fields(as.vector(tr))
  expect_true(length(tf$n) > 0 && all(tf$n == 11), info = "each transposon has 11 '!'-separated fields")
  if (!is.null(tf$m)) {
    expect_true(all(grepl("^[gi][01]MOCK_[0-9]{4}$", tf$m[, 1])), info = "gene field = g|i + copy 0|1 + locus tag")
    expect_true(all(tf$m[, 5] == uniq), info = "uniqID")
    expect_true(all(as.numeric(tf$m[, 6]) == tpn_tmpl$size[match(uniq, tpn_tmpl$uniqID)]), info = "size")
  }
  future <- vapply(seq_len(nrow(tr)), function(r) {
    m <- tpn_fields(tr[r, ])$m
    !is.null(m) && any(as.numeric(m[, 3]) > r - 1)
  }, logical(1))
  expect_false(any(future), info = "a transposon cannot carry a later generation than its record")
  # homologous repair restores only the jumped-from site and recombination replaces the recipient's segment,
  # so no host should ever hold the same transposon twice
  dups <- same_spot_duplicates(tr)
  expect_true(length(dups) == 0, info = paste(length(dups), "host records repeat a transposon, e.g.", dups[1]))

  g <- e$rec.gRecom
  expect_true(all(is.finite(g) & g >= 0))
  if (!is.null(recom_rate)) expect_true(all(abs(g - recom_rate) < 1e-12), info = "fixed homologous repair rate")

  # the analysis runs inside the generation loop; the saved one must cover every recorded generation
  if (file.exists(run$ana)) {
    a <- load_rda(run$ana)
    expect_true(all(c("gEnealogy", "rec.hostTraject", "tpn.dom", "tpn.gene", "g.Fit", "g.Rec", "g.Fit.ref") %in% ls(a)))
    expect_equal(unname(rowSums(a$gEnealogy)), rep(1, gmax + 1), tolerance = 1e-8)
    expect_equal(unname(rowSums(a$rec.hostTraject)), rep(1, gmax + 1), tolerance = 1e-8)
  }
  invisible(e)
}

expect_scenario <- function(name) {
  it <- integration()
  skip_if_not(isTRUE(it$ok[[name]]), paste("no template rows match scenario", name))
  k <- it$index[[name]]
  h <- it$tab$host[k]
  run <- run_sim(sb_main, k)
  fixed_rate <- if (identical(host_tmpl$hAReomH1[h], "fixed")) host_tmpl$homologousAutoRecom[h] else NULL
  expect_sim_outputs(run, it$tab$transposon[k], host_tmpl$cell[h], recom_rate = fixed_rate)
  invisible(run)
}

test_that("end to end: simulate.r's default scenario", {
  skip_integration()
  expect_scenario("default")
})

test_that("end to end: fixed rates, haploid host, recombination and homologous repair", {
  skip_integration()
  expect_scenario("fixed")
})

test_that("end to end: diploid host", {
  skip_integration()
  expect_scenario("diploid")
})

test_that("end to end: charlesworth jump rate and evolving copy rate", {
  skip_integration()
  expect_scenario("h1_rules")
})

test_that("end to end: genotoxic perturbation generations", {
  skip_integration()
  expect_scenario("genotoxic")
})

test_that("end to end: transposons that copy but never jump (jumpRate 0)", {
  skip_integration()
  expect_scenario("no_jumps")
})

test_that("end to end: copies restricted to the terminus direction", {
  skip_integration()
  expect_scenario("copy_dir")
})

test_that("end to end: copy direction in generations with copies but no jumps", {
  skip_integration()
  expect_scenario("copy_dir_copies_only")
})

test_that("end to end: copy direction in generations with jumps but no copies", {
  skip_integration()
  expect_scenario("copy_dir_jumps_only")
})

test_that("end to end: very low recombination rate (recom 0.0001)", {
  skip_integration()
  expect_scenario("recom_tiny")
})

test_that("end to end: recombination replaces the recipient's segment (repair off)", {
  skip_integration()
  run <- expect_scenario("recomb_only")     # includes the same-spot duplicate check
  skip_if(is.null(run) || run$status != 0 || !file.exists(run$sim), "run failed")
  expect_true(any(nzchar(load_rda(run$sim)$rec.host[-1, ])),
              info = "no recombination was recorded, so the duplicate check was not exercised")
})

test_that("end to end: homologous repair restores only the jumped-from site (recombination off)", {
  skip_integration()
  run <- expect_scenario("repair_only")     # includes the same-spot duplicate check
  skip_if(is.null(run) || run$status != 0 || !file.exists(run$sim), "run failed")
  m <- tpn_fields(as.vector(load_rda(run$sim)$rec.transposon[-1, ]))$m
  expect_true(!is.null(m) && any(as.numeric(m[, 3]) > 0),
              info = "no transposon moved, so repairs were not exercised")
})

test_that("end to end: offspring carry the transposons of the parent genomes recorded in familyTree", {
  skip_integration()
  run <- expect_scenario("inert")           # transposons never move and the host never recombines
  skip_if(is.null(run) || run$status != 0 || !file.exists(run$sim), "run failed")
  e <- load_rda(run$sim)
  bad <- character(0)
  for (r in 2:nrow(e$rec.offspring)) {
    for (j in seq_len(ncol(e$rec.offspring))) {
      k <- as.integer(strsplit(e$rec.offspring[r, j], ";", fixed = TRUE)[[1]])
      child <- chrom_keys(e$rec.transposon[r, j])
      for (slot in 1:2) {
        h <- (k[slot] - 1) %% HPOP + 1
        cp <- (k[slot] - 1) %/% HPOP
        if (!identical(child[[slot]], chrom_keys(e$rec.transposon[r - 1, h])[[cp + 1]])) {
          bad <- c(bad, sprintf("generation %d host %d copy %d vs parent genome %d (host %d, copy %d)",
                                r - 1, j, slot - 1, k[slot], h, cp))
        }
      }
    }
  }
  expect_true(length(bad) == 0, info = paste(length(bad), "mismatches, e.g.", bad[1]))
})

test_that("end to end: recombination copies donor transposons into haploid recipients", {
  skip_integration()
  run <- expect_scenario("hgt")             # inert transposons, haploid host recombining, repair off
  skip_if(is.null(run) || run$status != 0 || !file.exists(run$sim), "run failed")
  e <- load_rda(run$sim)
  gained <- 0
  foreign <- character(0)
  for (r in 2:nrow(e$rec.offspring)) {
    pool <- unique(unlist(lapply(e$rec.transposon[r - 1, ], function(cell) unlist(copy_keys(cell)))))
    for (j in seq_len(ncol(e$rec.offspring))) {
      k <- as.integer(strsplit(e$rec.offspring[r, j], ";", fixed = TRUE)[[1]])
      new <- setdiff(unlist(copy_keys(e$rec.transposon[r, j])), genome_keys(e$rec.transposon[r - 1, ], k[1]))
      gained <- gained + length(new)
      foreign <- c(foreign, setdiff(new, pool))
    }
  }
  expect_gt(gained, 0)                      # transposons the parent genome did not carry arrived from donors
  expect_true(length(foreign) == 0, info = paste("transposons found in no host of the previous generation, e.g.",
                                                 foreign[1]))
})

test_that("end to end: gene conversion copies the homolog's transposons in diploid hosts", {
  skip_integration()
  run <- expect_scenario("conversion")      # inert transposons, diploid host recombining
  skip_if(is.null(run) || run$status != 0 || !file.exists(run$sim), "run failed")
  e <- load_rda(run$sim)
  gained <- 0
  foreign <- character(0)
  for (r in 2:nrow(e$rec.offspring)) {
    for (j in seq_len(ncol(e$rec.offspring))) {
      k <- as.integer(strsplit(e$rec.offspring[r, j], ";", fixed = TRUE)[[1]])
      parents <- lapply(k, function(g) genome_keys(e$rec.transposon[r - 1, ], g))
      child <- copy_keys(e$rec.transposon[r, j])
      for (cp in 1:2) {
        new <- setdiff(child[[cp]], parents[[cp]])
        gained <- gained + length(new)
        foreign <- c(foreign, setdiff(new, parents[[3 - cp]]))   # must come from the homolog
      }
    }
  }
  expect_gt(gained, 0)
  expect_true(length(foreign) == 0, info = paste("transposons not taken from the homolog, e.g.", foreign[1]))
})

test_that("end to end: evolved rates are carried by transposons that move; unmoved ones keep the template rate", {
  skip_integration()
  it <- integration()
  skip_if_not(isTRUE(it$ok[["h1_rules"]]), "no template row with an evolving copy rate")
  k <- it$index[["h1_rules"]]
  run <- run_sim(sb_main, k)
  skip_if(run$status != 0 || !file.exists(run$sim), "the h1_rules run failed (see its own test)")
  m <- tpn_fields(as.vector(load_rda(run$sim)$rec.transposon[-1, ]))$m
  skip_if(is.null(m), "no transposons recorded")
  rate <- tpn_tmpl$copyRate[match(it$tab$transposon[k], tpn_tmpl$uniqID)]
  moved <- as.numeric(m[, 3]) > 0                       # generation stamp > 0: jumped or copied since the start
  expect_true(any(abs(as.numeric(m[moved, 9]) - rate) > 1e-9), info = "no moved transposon carries an evolved rate")
  expect_true(all(abs(as.numeric(m[!moved, 9]) - rate) < 1e-12), info = "an unmoved transposon's rate changed")
})

test_that("end to end: chromosome copies left without transposons (one transposon per host)", {
  skip_integration()
  run <- run_sim(sb_extra("sparse"), 1)
  expect_equal(run$status, 0, info = run$tail)
  expect_true(any(grepl("simulation completed", run$log, fixed = TRUE)), info = run$tail)
  if (run$status == 0 && file.exists(run$sim)) {
    dups <- same_spot_duplicates(load_rda(run$sim)$rec.transposon)
    expect_true(length(dups) == 0, info = paste(length(dups), "host records repeat a transposon, e.g.", dups[1]))
  }
})

test_that("end to end: a host population eliminated in generation 1 finishes cleanly", {
  skip_integration()
  run <- run_sim(sb_extra("lethal"), 1)
  expect_equal(run$status, 0, info = run$tail)
  expect_true(any(grepl("Host population eliminated", run$log, fixed = TRUE)), info = run$tail)
  # only generation 0 is recorded, so simulate.r runs no analysis; any analysis written must still be consistent
  if (file.exists(run$ana)) {
    a <- load_rda(run$ana)
    expect_equal(unname(rowSums(a$gEnealogy)), rep(1, nrow(a$gEnealogy)), tolerance = 1e-8)
  }
})

test_that("end to end: simulate.r reads the parameter file in its first argument, wherever it is", {
  skip_integration()
  it <- integration()
  skip_if_not(isTRUE(it$ok[["fixed"]]), "no template rows for the fixed scenario")
  k <- it$index[["fixed"]]
  h <- it$tab$host[k]
  alt_dir <- file.path(sb_main$root, "alt")
  dir.create(alt_dir, showWarnings = FALSE)
  alt <- utils::read.csv(file.path(sb_main$raw, "input.csv"), header = TRUE)
  alt$Value[alt$Type == "host organism constant population size"] <- "12"
  alt$Value[alt$Type == "host organism constant generation number"] <- "2"
  utils::write.csv(alt, file.path(alt_dir, "params.csv"), row.names = FALSE)
  run <- run_sim(sb_main, k, seed_index = 2L, input = "../alt/params.csv")   # templates still come from ../raw
  expect_sim_outputs(run, it$tab$transposon[k], host_tmpl$cell[h], gmax = 2, hpop = 12,
                     recom_rate = host_tmpl$homologousAutoRecom[h])
})
