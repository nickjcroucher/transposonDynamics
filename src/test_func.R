#!/usr/bin/env Rscript

library(testthat)

func_path <- Sys.getenv(
  "FUNC_R_PATH",
  unset = file.path(Sys.getenv("GITHUB_WORKSPACE", unset = getwd()), "src", "func.r")
)

if (!file.exists(func_path)) {
  stop("Cannot find func.r. Set FUNC_R_PATH to the script path before running these tests.")
}

source(normalizePath(func_path, mustWork = TRUE))

project_dir <- normalizePath(file.path(dirname(func_path), ".."), mustWork = TRUE)
src_dir <- dirname(func_path)

with_dir <- function(path, code) {
  old <- getwd()
  on.exit(setwd(old), add = TRUE)
  setwd(path)
  force(code)
}

toy_gene_df <- function() {
  data.frame(
    locus_tag = c("A", "B", "C"),
    start = c(1, 101, 201),
    end = c(100, 200, 300),
    product = c("prodA", "prodB", "prodC"),
    length = c(100, 100, 100),
    interLength = c(0, 10, 20),
    essential = c(FALSE, FALSE, FALSE),
    recombination = c(FALSE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
}

test_that("gffClean expands GFF attributes into columns", {
  gff <- data.frame(
    V1 = c("chr1", "chr1"),
    V2 = c("src", "src"),
    V3 = c("region", "CDS"),
    V4 = c(1, 5),
    V5 = c(100, 9),
    V6 = c(".", "."),
    V7 = c("+", "+"),
    V8 = c(".", "0"),
    V9 = c(
      "ID=region1;Name=chromosome",
      "ID=cds1;locus_tag=geneA;product=proteinA"
    ),
    stringsAsFactors = FALSE
  )

  clean <- gffClean(gff)

  expect_named(
    clean,
    c(
      "seqname", "source", "type", "start", "end", "score", "strand",
      "phase", "ID", "Name", "locus_tag", "product"
    )
  )
  expect_equal(clean$locus_tag[2], "geneA")
  expect_equal(clean$product[2], "proteinA")
  expect_true(is.na(clean$product[1]))
})

test_that("rNumVec delegates to normal, uniform, and poisson generators", {
  set.seed(123)
  actual_normal <- rNumVec("normal", L = "3", p1 = "10", p2 = "2")
  set.seed(123)
  expected_normal <- rnorm(n = 3, mean = 10, sd = 2)
  expect_equal(actual_normal, expected_normal)

  set.seed(123)
  actual_uniform <- rNumVec("uniform", L = 3, p1 = 2, p2 = 5)
  set.seed(123)
  expected_uniform <- runif(n = 3, min = 0, max = 7)
  expect_equal(actual_uniform, expected_uniform)

  set.seed(123)
  actual_poisson <- rNumVec("poisson", L = 3, p1 = .3)
  set.seed(123)
  expected_poisson <- 1/rpois(n = 3, lambda = .3)
  expect_equal(actual_poisson, expected_poisson)

  set.seed(123)
  actual_poisson <- rNumVec("negbin", L = 3, p1 = .3, p2 = .4)
  set.seed(123)
  expected_poisson <- 1/rnbinom(n = 3, prob = .3, size = .4)
  expect_equal(actual_poisson, expected_poisson)

  expect_error(
    rNumVec("negative-binomial", L = 1, p1 = 1, p2 = 1),
    "Only three function options allowed"
  )
})

test_that("inParams loads parameters, gene table, and initial transposon titres", {
  skip_if_not(file.exists(file.path(project_dir, "raw", "input.csv")))
  skip_if_not(file.exists(file.path(project_dir, "data", "GCA_001457635.gff")))

  set.seed(42)
  out <- with_dir(src_dir, inParams("../raw/input.csv"))
  pop_size <- as.numeric(
    out$params$Value[out$params$Type == "host organism constant population size"]
  )

  expect_named(out, c("params", "gene", "transposon.titre"))
  expect_true(all(c(
    "locus_tag", "start", "end", "product", "length", "interLength",
    "essential", "recombination"
  ) %in% names(out$gene)))
  expect_gt(nrow(out$params), 0)
  expect_gt(nrow(out$gene), 0)
  expect_length(out$transposon.titre, pop_size)
  expect_true(all(out$gene$length > 0))
  expect_true(all(out$gene$essential == FALSE))
  expect_true(all(out$gene$recombination == FALSE))
})

test_that("ini.host creates compact host genomes and validates variation counts", {
  gene_df <- toy_gene_df()

  set.seed(10)
  hosts <- ini.host(host.var = 6, gene.df = gene_df, gene.var = "1;2;3")
  host_chars <- do.call(rbind, strsplit(hosts, split = ""))

  expect_length(hosts, 6)
  expect_true(all(nchar(hosts) == nrow(gene_df)))
  expect_true(all(host_chars[, 1] == gVar[1]))
  expect_true(all(host_chars[, 2] %in% gVar[1:2]))
  expect_true(all(host_chars[, 3] %in% gVar[1:3]))

  one_gene_hosts <- ini.host(
    host.var = 3,
    gene.df = gene_df[1, , drop = FALSE],
    gene.var = "1"
  )
  expect_equal(one_gene_hosts, rep("a", 3))

  expect_warning(
    too_many <- ini.host(host.var = 2, gene.df = gene_df, gene.var = "1;1;1;1"),
    "Too many notations"
  )
  expect_equal(too_many, rep("aaa", 2))

  expect_error(
    ini.host(host.var = 2, gene.df = gene_df, gene.var = "1;2"),
    "Either one variation or specify variations for each gene"
  )
})

test_that("ini.transposon creates encoded transposon records", {
  set.seed(11)
  out <- ini.transposon("100;200")
  parts <- strsplit(out$ini, split = "!", fixed = TRUE)

  expect_s3_class(out, "data.frame")
  expect_named(out, c("ini", "uniqID", "size"))
  expect_equal(out$size, c(100, 200))
  expect_match(out$uniqID, "^[A-Z]{7}$")
  expect_equal(vapply(parts, `[[`, character(1), 1), rep("", 2))
  expect_equal(vapply(parts, `[[`, character(1), 2), rep("0", 2))
  expect_equal(vapply(parts, `[[`, character(1), 3), rep("T", 2))
  expect_equal(vapply(parts, `[[`, character(1), 4), out$uniqID)
})

test_that("reZero rescales values with the provided bounds", {
  expect_equal(reZero(5, new0 = 0, new1 = 10), 0.5)
  expect_equal(reZero(c(2, 4, 6), new0 = 2, new1 = 6), c(0, 0.5, 1))
})

test_that("reGeneDF shifts gene coordinates after transposon insertion", {
  gene_df <- toy_gene_df()

  gene_insert <- reGeneDF("gA!5!55!T!id1", tPn.size = 10, gene.df = gene_df)
  expect_equal(gene_insert$length, c(110, 100, 100))
  expect_equal(gene_insert$start, c(1, 111, 211))
  expect_equal(gene_insert$end, c(110, 210, 310))

  intergenic_insert <- reGeneDF("iB!5!55!T!id2", tPn.size = 5, gene.df = gene_df)
  expect_equal(intergenic_insert$interLength, c(0, 15, 20))
  expect_equal(intergenic_insert$start, c(1, 106, 206))
  expect_equal(intergenic_insert$end, c(100, 205, 305))
})

test_that("host.reproduce samples offspring and excludes essential-gene insertions", {
  gene_df <- toy_gene_df()
  res_pool <- data.frame(
    host = c("aaa", "bbb", "ccc"),
    transposon = c("gA!1!5!T!id1", "gB!2!5!T!id2", "gC!3!5!T!id3"),
    stringsAsFactors = FALSE
  )

  set.seed(1)
  out <- host.reproduce(res.pool = res_pool, gene.df = gene_df, transposon.size = 0)
  set.seed(1)
  expected_idx <- sample(
    seq_len(nrow(res_pool)),
    nrow(res_pool),
    replace = TRUE,
    prob = rep(1 / nrow(res_pool), nrow(res_pool))
  )

  expect_equal(out$familyTree, expected_idx)
  expect_equal(out$host, res_pool$host[expected_idx])
  expect_equal(out$transposon, res_pool$transposon[expected_idx])

  gene_df$essential[2] <- TRUE
  set.seed(1)
  filtered <- host.reproduce(res.pool = res_pool, gene.df = gene_df, transposon.size = 0)
  expect_false(any(filtered$familyTree == 2))
})

test_that("tPn.io converts between fields and encoded chains", {
  encoded <- tPn.io(gene = "gA", location = 12, generation = 2, valid = "T", uniqID = "id1")

  expect_equal(encoded, "gA!12!2!T!id1")
  expect_equal(
    tPn.io(chain = encoded, encrypt = FALSE),
    c("gA", "12", "2", "T", "id1")
  )
})

test_that("tPn.x invalidates the earlier overlapping transposon", {
  overlap <- tPn.x(
    tPn1 = "gA!5!5!T!id1",
    tPn1.len = 3,
    tPn2 = "gB!7!5!T!id2",
    tPn2.len = 1
  )
  expect_equal(overlap, "gA!5!5!F!id1")

  no_overlap <- tPn.x(
    tPn1 = "gA!5!5!T!id1",
    tPn1.len = 3,
    tPn2 = "gB!9!5!T!id2",
    tPn2.len = 1
  )
  expect_equal(no_overlap, "gA!5!5!T!id1")

  already_invalid <- tPn.x(
    tPn1 = "gA!5!5!F!id1",
    tPn1.len = 3,
    tPn2 = "gB!7!5!T!id2",
    tPn2.len = 1
  )
  expect_equal(already_invalid, c("gA!5!5!F!id1", "gB!7!5!T!id2"))
})

test_that("tPn.r sets the fourth field to T for one record", {
  expect_equal(
    tPn.r("gA!12!0!T!id001"),
    "gA!12!0!T!id001"
  )
})

test_that("tPn.r sets the fourth field to T for multiple records", {
  expect_equal(
    tPn.r("gA!12!0!T!id001;iB!34!23!F!id002"),
    "gA!12!0!T!id001;iB!34!23!T!id002"
  )
})

test_that("tPn.r preserves record count and the first three fields", {
  input <- "gA!12!0!T!id001;iB!34!23!F!id002"
  output <- tPn.r(input)

  in_parts <- strsplit(strsplit(input, ";", fixed = TRUE)[[1]], "!", fixed = TRUE)
  out_parts <- strsplit(strsplit(output, ";", fixed = TRUE)[[1]], "!", fixed = TRUE)

  expect_length(out_parts, length(in_parts))
  expect_equal(lapply(out_parts, `[`, 1:3), lapply(in_parts, `[`, 1:3))
  expect_equal(vapply(out_parts, `[`, character(1), 4), rep("T", length(out_parts)))
})

test_that("tPn.jump leaves tags unchanged when no jump happens", {
  gene_df <- toy_gene_df()
  tag <- "gA!1!5!T!id1"

  out <- tPn.jump(
    tPn.tag = tag,
    tPn.prob = 0.9,
    tPn.size = 5,
    gene.df = gene_df,
    jumProb = 0.1
  )

  expect_equal(out$tPn.tag, tag)
  expect_equal(out$gene, gene_df)
})

test_that("tPn.jump returns an adjusted gene table after a jump", {
  gene_df <- toy_gene_df()

  out <- tPn.jump(
    tPn.tag = "gA!1!5!T!id1",
    tPn.prob = 0.5,
    tPn.size = 5,
    gene.df = gene_df,
    tPn.move = 0.5,
    jumProb = 1,
    generation = 8
  )
  fields <- tPn.io(chain = out$tPn.tag, encrypt = FALSE)

  expect_length(fields, 5)
  expect_true(substr(fields[1], 1, 1) %in% c("g", "i"))
  expect_true(substring(fields[1], 2) %in% gene_df$locus_tag)
  expect_true(as.numeric(fields[2]) > 0)
  expect_true(as.numeric(fields[3]) >= 0)
  expect_equal(fields[4], "T")
  expect_equal(fields[5], "id1")
  expect_false("cumsum" %in% names(out))
  expect_equal(nrow(out), nrow(gene_df))
  expect_true(any(out$end > gene_df$end))
})

test_that("gene.recom replaces one allele and transfers matching transposons", {
  out <- gene.recom(
    h1.G = "aaa",
    h1.t = "gB!5!2!T!h1B;gC!1!3!T!h1C",
    h2.G = "BCD",
    h2.t = "gB!7!6!T!h2B;gA!2!5!T!h2A",
    g2to1 = 2,
    locusTags = c("A", "B", "C")
  )

  expect_equal(unname(out["genome"]), "aCa")
  expect_equal(unname(out["transposon"]), "gC!1!3!T!h1C;gB!7!6!T!h2B")
})

test_that("g.Recom leaves rows untouched with zero recombination and changes recipients otherwise", {
  gene_df <- toy_gene_df()
  res_pool <- data.frame(
    host = c("aaa", "BBB", "CCC"),
    transposon = c("gA!1!2!T!a1", "gB!2!4!T!b1", "gC!3!6!T!c1"),
    stringsAsFactors = FALSE
  )

  expect_equal(suppressWarnings(g.Recom(res_pool, gene_df, c(0, 0, 0))), res_pool)

  set.seed(1)
  out <- g.Recom(res_pool, gene_df, c(1, 0, 0))

  expect_equal(out$host[2:3], res_pool$host[2:3])
  expect_equal(out$transposon[2:3], res_pool$transposon[2:3])
  expect_equal(nchar(out$host), rep(nrow(gene_df), nrow(out)))
  expect_equal(
    sum(strsplit(out$host[1], split = "")[[1]] != strsplit(res_pool$host[1], split = "")[[1]]),
    1
  )
})
