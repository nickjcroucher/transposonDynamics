# Unit tests for func.r.
#
# func.r is source()'d fresh (into its own environment, `fenv`) from an
# isolated sandbox copy for every test file run -- see
# tests/testthat/helper-fixtures.R for how that sandbox is built and why.
#
# Most tests below use a small hand-built gene table (make_gene_df()) rather
# than a real genome, because inParams()/gffClean() need an external genome
# GFF file that is not part of this repo's four supplied CSVs (see the
# dedicated inParams() tests near the end of this file, which skip
# gracefully when no GFF is present).
#
# A few tests are deliberately written to assert func.r's CURRENT behavior
# around some fragile-looking code, rather than "corrected" behavior. Those
# are commented as such; they document/highlight the issue without fixing
# func.r itself, per the review request.

library(testthat)

sandbox <- build_sandbox()
fenv <- source_func_r(sandbox)

## ---- gffClean() ------------------------------------------------------

test_that("gffClean parses ';'-separated GFF attributes into their own columns", {
  raw_gff <- data.frame(
    V1 = c("chr1", "chr1", "chr1"),
    V2 = c("src", "src", "src"),
    V3 = c("region", "gene", "gene"),
    V4 = c(1, 1, 501),
    V5 = c(1000, 400, 900),
    V6 = c(".", ".", "."),
    V7 = c("+", "+", "+"),
    V8 = c(".", ".", "."),
    V9 = c(
      "ID=chr1",
      "ID=gene1;locus_tag=SPD_0001;product=hypothetical protein A",
      "ID=gene2;locus_tag=SPD_0002;product=protein B"
    ),
    stringsAsFactors = FALSE
  )
  out <- fenv$gffClean(raw_gff)
  expect_equal(nrow(out), 3)
  expect_true(all(c("locus_tag", "product", "start", "end", "type") %in% colnames(out)))
  expect_equal(out$locus_tag[2], "SPD_0001")
  expect_equal(out$product[3], "protein B")
  expect_true(is.na(out$locus_tag[1])) # the "region" row has no locus_tag attribute
})

test_that("gffClean silently drops attribute text after a second '=' (documents func.r:18-24 behavior)", {
  tricky <- data.frame(
    V1 = "chr1", V2 = "src", V3 = "gene", V4 = 1, V5 = 100,
    V6 = ".", V7 = "+", V8 = ".",
    V9 = "ID=gene9;locus_tag=SPD_0009;product=a=b",
    stringsAsFactors = FALSE
  )
  out <- fenv$gffClean(tricky)
  # gffClean splits each "key=value" piece on EVERY '=' and keeps only the
  # second piece (func.r line 23-24), so a value that itself contains '='
  # loses everything after the second '=' -- here "a=b" becomes just "a".
  # This is flagged as a warning in the review, not corrected here.
  expect_equal(out$product, "a")
})

## ---- rNumVec() ---------------------------------------------------------

test_that("rNumVec's four distributions all return a vector of length L", {
  set.seed(1)
  expect_length(fenv$rNumVec("normal", L = 5, p1 = 0, p2 = 1), 5)
  expect_length(fenv$rNumVec("uniform", L = 5, p1 = 0, p2 = 1), 5)
  expect_length(fenv$rNumVec("poisson", L = 5, p1 = 2), 5)
  expect_length(fenv$rNumVec("negbin", L = 5, p1 = 0.5, p2 = 2), 5)
})

test_that("rNumVec's poisson/negbin branches always return positive values", {
  set.seed(2)
  expect_true(all(fenv$rNumVec("poisson", L = 200, p1 = 3) > 0))
  expect_true(all(fenv$rNumVec("negbin", L = 200, p1 = 0.3, p2 = 4) > 0))
})

test_that("rNumVec errors on an unknown distribution, and the message omits 'negbin' (func.r:37)", {
  expect_error(fenv$rNumVec("bogus", L = 1), "normal, uniform, poisson")
  err <- tryCatch(fenv$rNumVec("bogus", L = 1), error = function(e) conditionMessage(e))
  # negbin (func.r:36) is a valid 4th option that the error text never mentions
  expect_false(grepl("negbin", err))
})

test_that("inParams' transposon-titre formula is only sensible for poisson/negbin, not normal/uniform (func.r:58)", {
  # inParams() computes: round(1 / rNumVec(...), 0)  [func.r line 58]
  # rNumVec's poisson/negbin branches already return 1/(x+1) internally
  # (func.r lines 35-36), so for those two the outer 1/(...) cancels back
  # out to a sensible count. normal/uniform have no internal inversion, so
  # the SAME outer 1/(...) is left un-cancelled -- a typical mean like 5
  # silently collapses to 0 for every host. All four are equally valid,
  # documented choices in input.csv's Note column.
  set.seed(21)
  poisson_titre <- round(1 / fenv$rNumVec("poisson", L = 1000, p1 = 5, p2 = 0), 0)
  normal_titre <- round(1 / fenv$rNumVec("normal", L = 1000, p1 = 5, p2 = 0), 0)
  expect_true(mean(poisson_titre) > 1) # sensible: averages close to mean+1 = 6
  expect_true(mean(normal_titre) == 0) # collapses to 0 for every host
})

## ---- reZero() ------------------------------------------------------------

test_that("reZero linearly rescales a value (or vector) between two bounds", {
  expect_equal(fenv$reZero(5, new0 = 0, new1 = 10), 0.5)
  expect_equal(fenv$reZero(c(0, 5, 10), new0 = 0, new1 = 10), c(0, 0.5, 1))
})

test_that("reZero has no guard against new0 == new1 (func.r:100) and returns Inf/NaN instead of erroring", {
  result <- fenv$reZero(5, new0 = 3, new1 = 3)
  expect_true(is.nan(result) || is.infinite(result))
})

## ---- ini.host() ----------------------------------------------------------

test_that("ini.host builds one genome string per host, one symbol per gene", {
  gdf <- make_gene_df()
  set.seed(1)
  hosts <- fenv$ini.host(host.var = 5, gene.df = gdf, gene.var = "10")
  expect_length(hosts, 5)
  expect_true(all(nchar(hosts) == nrow(gdf)))
})

test_that("ini.host warns and truncates when given more per-gene variation notations than genes", {
  gdf <- make_gene_df()
  set.seed(2)
  expect_warning(
    hosts <- fenv$ini.host(host.var = 4, gene.df = gdf, gene.var = "5;10;3;99"),
    "Too many notations"
  )
  expect_length(hosts, 4)
  expect_true(all(nchar(hosts) == nrow(gdf)))
})

test_that("ini.host errors when given more than one but too few variation notations", {
  gdf <- make_gene_df()
  expect_error(
    fenv$ini.host(host.var = 4, gene.df = gdf, gene.var = "5;10"),
    "Either one variation"
  )
})

## ---- ini.transposon() -----------------------------------------------------

test_that("ini.transposon builds one record per requested size, using the template's column count", {
  set.seed(3)
  scenario_row <- data.frame(jumpRate = 0.1, copyRate = 0.2)
  out <- fenv$ini.transposon(tPn.size = "1000;1500", scenario = scenario_row, template = fenv$tPn.0)
  expect_equal(nrow(out), 2)
  expect_true(all(c("ini", "uniqID", "size") %in% colnames(out)))
  expect_equal(out$size, c(1000, 1500))
  expect_equal(unname(lengths(strsplit(out$ini, "!"))), rep(ncol(fenv$tPn.0), 2))
})

## ---- reGeneDF() ------------------------------------------------------------

test_that("reGeneDF grows gene length for a coding-region ('g') insertion and shifts downstream start/end", {
  gdf <- make_gene_df()
  out <- fenv$reGeneDF(tPn = "gSPD_0001!50!0!TRUE!AAAAAAA!1000!0.1!0.1", tPn.size = 1000, gene.df = gdf)
  expect_equal(out$length[1], gdf$length[1] + 1000)
  expect_equal(out$start[2], gdf$start[2] + 1000)
  expect_equal(out$end[1], gdf$end[1] + 1000)
  expect_equal(out$end[3], gdf$end[3] + 1000)
})

test_that("reGeneDF grows intergenic length for an intergenic ('i') insertion", {
  gdf <- make_gene_df()
  out <- fenv$reGeneDF(tPn = "iSPD_0002!20!0!TRUE!BBBBBBB!1000!0.1!0.1", tPn.size = 1000, gene.df = gdf)
  expect_equal(out$interLength[2], gdf$interLength[2] + 1000)
  expect_equal(out$start[2], gdf$start[2] + 1000)
  expect_equal(out$end[2], gdf$end[2] + 1000)
})

test_that("reGeneDF is a no-op when there is no transposon", {
  gdf <- make_gene_df()
  out <- fenv$reGeneDF(tPn = "", tPn.size = 1000, gene.df = gdf)
  expect_equal(out, gdf)
})

## ---- host.reproduce() -------------------------------------------------------

test_that("host.reproduce gives zero reproduction probability to a host with a transposon in an essential gene", {
  gdf <- make_gene_df() # SPD_0001 (gene 1) is essential
  res_pool <- data.frame(
    host = c("AAA", "BBB", "CCC", "DDD"),
    transposon = c("gSPD_0001!10!0!TRUE!XXXXXXX!1000!0.1!0.1", "", "", ""),
    stringsAsFactors = FALSE
  )
  out <- fenv$host.reproduce(res.pool = res_pool, gene.df = gdf, transposon.size = 1000, fitness.advantage = 10)
  # a probability of exactly zero makes this deterministic, not just likely
  expect_false(1 %in% out$familyTree)
  expect_equal(nrow(out), 4)
})

## ---- tPn.io() ------------------------------------------------------------

test_that("tPn.io converts a value vector to a '!'-joined string and back to a named vector", {
  vals <- c("gSPD_0001", "10", "0", "TRUE", "ABCDEFG", "1000", "0.1", "0.1")
  s <- fenv$tPn.io(vals)
  expect_type(s, "character")
  expect_length(strsplit(s, "!")[[1]], 8)

  named <- fenv$tPn.io(s)
  expect_equal(unname(named["gene"]), "gSPD_0001")
  expect_equal(unname(named["size"]), "1000")
})

test_that("tPn.io converts a ';'-joined multi-transposon string into a data.frame", {
  s1 <- fenv$tPn.io(c("gSPD_0001", "10", "0", "TRUE", "ABCDEFG", "1000", "0.1", "0.1"))
  s2 <- fenv$tPn.io(c("iSPD_0002", "5", "0", "TRUE", "HIJKLMN", "1000", "0.2", "0.2"))
  df <- fenv$tPn.io(paste(s1, s2, sep = ";"))
  expect_equal(nrow(df), 2)
  expect_equal(ncol(df), 8)
})

test_that("tPn.io errors on a value vector of the wrong length", {
  expect_error(fenv$tPn.io(c("only", "two")))
})

## ---- tPn.x() -------------------------------------------------------------

test_that("tPn.x deactivates tPn1 when tPn2 overlaps it at the same locus", {
  t1 <- "gSPD_0001!2!0!TRUE!ABCDEFG!3!0.1!0.1"
  t2 <- "gSPD_0001!1!0!TRUE!ZZZZZZZ!6!0.1!0.1"
  out <- fenv$tPn.x(t1, t2)
  expect_equal(strsplit(out, "!")[[1]][4], "FALSE")
})

test_that("tPn.x leaves tPn1 unchanged when the two loci differ", {
  t1 <- "gSPD_0001!2!0!TRUE!ABCDEFG!3!0.1!0.1"
  t2 <- "gSPD_0002!1!0!TRUE!ZZZZZZZ!6!0.1!0.1"
  out <- fenv$tPn.x(t1, t2)
  expect_equal(out, t1)
})

## ---- tPn.r() -------------------------------------------------------------

test_that("tPn.r revives all transposons by forcing the 4th (valid) field to TRUE", {
  s1 <- "gSPD_0001!10!0!FALSE!ABCDEFG!1000!0.1!0.1"
  s2 <- "iSPD_0002!5!0!FALSE!HIJKLMN!1000!0.2!0.2"
  out <- fenv$tPn.r(paste(s1, s2, sep = ";"))
  expect_length(out, 2)
  fields <- strsplit(out, "!")
  expect_true(all(vapply(fields, `[`, character(1), 4) == "TRUE"))
})

## ---- tPn.reloc() ---------------------------------------------------------

test_that("tPn.reloc grows total genome length by exactly the transposon size, wherever it lands", {
  gdf <- make_gene_df()
  total_before <- sum(gdf$length) + sum(gdf$interLength)
  set.seed(7)
  out <- fenv$tPn.reloc(tPn = "gSPD_0001!1!0!TRUE!AAAAAAA!1000!0.1!0.1", gen = 1, gene.df = gdf)
  total_after <- sum(out$length) + sum(out$interLength)
  expect_equal(total_after - total_before, 1000)
  expect_equal(nrow(out), nrow(gdf))
  expect_equal(ncol(out), ncol(gdf)) # the transient $cumsum column must be removed before returning
  expect_true(grepl(";", colnames(out)[1])) # new transposon tag prepended to the rolling history
})

## ---- h1.mod() ------------------------------------------------------------

test_that("h1.mod's default ('fixed rate') branch just repeats vAl for both outputs", {
  out <- fenv$h1.mod(0.2, numTpn = 3, H1 = "fixed rate", gToxic = 1)
  expect_equal(out, c(0.2, 0.2))
})

test_that("h1.mod's 'charlesworth' branch divides only the first (jump) value by numTpn", {
  out <- fenv$h1.mod(0.2, numTpn = 4, H1 = "charlesworth", gToxic = 1)
  expect_equal(out, c(0.05, 0.2))
})

test_that("h1.mod's 'evolving' branch returns two equal, randomly-perturbed non-negative values", {
  set.seed(42)
  out <- fenv$h1.mod(0.2, numTpn = 3, H1 = "evolving rate", gToxic = 1)
  expect_length(out, 2)
  expect_equal(out[1], out[2])
  expect_true(out[1] >= 0)
})

test_that("h1.mod applies gToxic only to the first (jump) element", {
  out <- fenv$h1.mod(0.2, numTpn = 3, H1 = "fixed rate", gToxic = 0.5)
  expect_equal(out, c(0.1, 0.2))
})

## ---- sCene.mod() ---------------------------------------------------------

test_that("sCene.mod reorders (jump, copy) pairs, and gToxic = FALSE bypasses the genotoxic branch", {
  sCene_row <- c(NA, "fixed", NA, "fixed") # only positions 2 and 4 (H1 strings) are read
  out <- fenv$sCene.mod(pRobs = c(0.3, 0.4), tPn.bg = "a;b;c", sCene = sCene_row, pAram = c(0.1, 0.1, 1), gToxic = FALSE)
  expect_equal(out, c(0.3, 0.4, 0.3, 0.4))
})

## ---- tPn.get() -----------------------------------------------------------

test_that("tPn.get splits the transposon-history colname into transposons vs. the current locus", {
  gene_df <- data.frame(x = 1)
  colnames(gene_df)[1] <- "tag1;tag2;currentLocusTag"
  expect_equal(fenv$tPn.get(gene_df, transposon = TRUE), "tag1;tag2")
  expect_equal(fenv$tPn.get(gene_df, transposon = FALSE), "currentLocusTag")
})

test_that("tPn.get returns an empty transposon string when none have been recorded yet", {
  gene_df <- data.frame(x = 1)
  colnames(gene_df)[1] <- "onlyLocusTag"
  expect_equal(fenv$tPn.get(gene_df, transposon = TRUE), "")
  expect_equal(fenv$tPn.get(gene_df, transposon = FALSE), "onlyLocusTag")
})

## ---- gene.recom() ---------------------------------------------------------

test_that("gene.recom transfers one gene and its transposon tag from host 2 into host 1", {
  out <- fenv$gene.recom(
    h1.G = "ABC", h1.t = "gSPD_0002!5!0!TRUE!AAAAAAA!1000!0.1!0.1",
    h2.G = "XYZ", h2.t = "gSPD_0002!7!0!TRUE!BBBBBBB!1000!0.2!0.2",
    g2to1 = 2, locusTags = c("SPD_0001", "SPD_0002", "SPD_0003")
  )
  expect_equal(out[["genome"]], "AYC")
  # h1's own SPD_0002 transposon is dropped and h2's is transferred in; note
  # gene.recom leaves a leading ';' artifact here -- its caller, g.Recom,
  # strips leading/trailing ';' itself (see func.r's g.Recom, last lines)
  expect_true(grepl("BBBBBBB", out[["transposon"]]))
  expect_false(grepl("AAAAAAA", out[["transposon"]]))
})

test_that("gene.recom leaves the transposon string untouched when neither host has a tag at that locus", {
  out <- fenv$gene.recom(
    h1.G = "ABC", h1.t = "gSPD_0001!5!0!TRUE!AAAAAAA!1000!0.1!0.1",
    h2.G = "XYZ", h2.t = "gSPD_0003!7!0!TRUE!BBBBBBB!1000!0.2!0.2",
    g2to1 = 2, locusTags = c("SPD_0001", "SPD_0002", "SPD_0003")
  )
  expect_equal(out[["genome"]], "AYC")
  expect_equal(out[["transposon"]], "gSPD_0001!5!0!TRUE!AAAAAAA!1000!0.1!0.1")
})

## ---- g.Recom() -----------------------------------------------------------

test_that("g.Recom excludes a host from recombining when its transposon hits a recombination-mechanism gene", {
  gdf <- make_gene_df() # SPD_0002 (gene 2) is flagged for recombination mechanism
  res_pool <- data.frame(
    host = c("AAAA", "TTTT", "CCCC"),
    transposon = c("gSPD_0002!5!0!TRUE!QQQQQQQ!1000!0.1!0.1", "", ""),
    stringsAsFactors = FALSE
  )
  set.seed(11)
  out <- fenv$g.Recom(res.pool = res_pool, gene.df = gdf, recomRate = 1, hypothesis = "switch")
  expect_equal(nrow(out), 3)
  expect_equal(ncol(out), ncol(res_pool))
  expect_equal(out$host[1], res_pool$host[1])
  expect_equal(out$transposon[1], res_pool$transposon[1])
})

test_that("g.Recom's per-host gene sampling excludes gene index == host row number, not a gene property (func.r's g.Recom)", {
  # sample(c(1:nrow(gene.df))[-i], numGenes[i], ...) uses the RES.POOL row
  # index i to exclude a GENE index -- two unrelated index spaces reused
  # under the same variable name. With 3 hosts and 3 genes this never
  # errors (there's always exactly nrow(gene.df)-1 genes left to sample),
  # but it means host i can structurally never recombine at gene i, for no
  # biological reason. This test just documents that it runs without error
  # for a same-sized pool; it does not assert the (questionable) exclusion
  # is intentional.
  gdf <- make_gene_df()
  res_pool <- data.frame(
    host = c("AAAA", "TTTT", "CCCC"),
    transposon = c("", "", ""),
    stringsAsFactors = FALSE
  )
  set.seed(12)
  expect_no_error(fenv$g.Recom(res.pool = res_pool, gene.df = gdf, recomRate = 1, hypothesis = "switch"))
})

## ---- inParams() (needs an external genome GFF not supplied with this repo) ----

test_that("inParams fails clearly when the referenced genome GFF is unavailable (documents an external dependency)", {
  sb <- build_sandbox()
  gffs <- list.files(file.path(sb, "data"), pattern = "\\.gff$")
  skip_if(length(gffs) > 0, "A real genome GFF is present in this repo; the missing-file path isn't reachable here.")
  fenv2 <- source_func_r(sb)
  old <- setwd(file.path(sb, "code"))
  on.exit(setwd(old), add = TRUE)
  expect_error(fenv2$inParams("../raw/input.csv"))
})

test_that("inParams succeeds end-to-end when a real genome GFF is available", {
  sb <- build_sandbox()
  gffs <- list.files(file.path(sb, "data"), pattern = "\\.gff$")
  skip_if_not(length(gffs) > 0, "No genome GFF found in data/ next to raw/ -- add one to exercise this end-to-end path.")
  fenv2 <- source_func_r(sb)
  old <- setwd(file.path(sb, "code"))
  on.exit(setwd(old), add = TRUE)
  result <- fenv2$inParams("../raw/input.csv")
  expect_true(all(c("params", "gene", "transposon.titre", "genome") %in% names(result)))
  expect_gt(nrow(result$gene), 0)
})
