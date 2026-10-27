#!/bin/env Rscript
# author: ph-u
# script: src.r
# desc: paths
# in: source("src.r")
# out: NA
# arg: 0
# date: 20261005

source("p_path.r")
suppressWarnings(library(readMDTable))
library(jsonlite)

src.url = "https://raw.githubusercontent.com/ph-u/snippet_ph-u/refs/heads/master/"
source(paste0(src.url, "colour.r"))
source(paste0(src.url, "eta.r"))

##### import and flatten JSON files #####
na_chr <- function(x) if (is.null(x)) NA_character_ else x

flat <- function(n, parent = NA_character_, depth = 0) {
  row <- data.frame(
    id = n$id, parent_id = parent, depth = depth,
    rank = na_chr(n$rank), name = na_chr(n$name),
    proteins = n$proteins, species = n$species)
  kids <- lapply(n$children, flat, n$id, depth + 1)
  do.call(rbind, c(list(row), kids))
}

##### taxa tree taxid #####
f.ref = read.table(f.ls[grep("table", f.ls)], sep = "\t", header = TRUE)

f.claude = as.data.frame(suppressMessages(read_md_table(f.ls[grep(".md$", f.ls)])))
colnames(f.claude) = read.table(text = sub("[.]","@",colnames(f.claude)), sep = "@")[,2]

system(paste("head -n 128", f.ls[grep("n.tsv", f.ls)], "| tail -n +7 > ../data/CL0219_gpt56.tsv"))
source("p_path.r")
f.gpt56 = read.table(f.ls[grep("gpt", f.ls)], sep = "\t", header = TRUE, comment.char = "", quote = "")
system("rm ../data/CL0219_gpt56.tsv")

f.0 = merge(f.claude[,c(1,3)], f.gpt56[,c(1,3)], all = TRUE, by = colnames(f.claude)[1])
