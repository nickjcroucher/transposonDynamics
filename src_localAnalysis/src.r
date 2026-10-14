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

f.ls = f.all[grep("CL0219", f.all)]
f.ref = read.table(f.ls[grep("table", f.ls)], sep = "\t", header = TRUE)

f.claude = as.data.frame(suppressMessages(read_md_table(f.ls[grep(".md$", f.ls)])))
colnames(f.claude) = read.table(text = sub("[.]","@",colnames(f.claude)), sep = "@")[,2]

f.gpt56 = read.table(f.ls[grep("gpt", f.ls)], sep = "\t", header = TRUE, comment.char = "", quote = "")
system("rm ../data/CL0219_gpt56.tsv")

f.0 = merge(f.claude[,c(1,3)], f.gpt56[,c(1,3)], all = TRUE, by = colnames(f.claude)[1])
