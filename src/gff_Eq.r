#!/bin/env Rscript
# author: ph-u
# script: gff_Eq.r
# desc: equate gff gene notations
# in: Rscript gff_Eq.r
# out: data/GCA_000014365.2.csv
# arg: 0
# date: 20260703

source("func.r")
a = gffClean(read.table("../data/GCA_000014365.gff", sep = "\t", header = F, quote = ""))
a0 = a[a$type=="CDS",]
a0$TIGR4 = NA
a0.x = grep(" TIGR4 ", a0$Note)
a0$TIGR4[a0.x] = substr(read.table(text = sub(" TIGR4 %3D ","@" , a0$Note[a0.x]), sep = "@")[,2], 1, 6)
write.csv(a0, "../data/GCA_000014365.csv", row.names = F, quote = F)

# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000014365.2/
# https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=1313
