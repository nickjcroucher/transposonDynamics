#!/bin/env Rscript
# author: ph-u
# script: 20260825_cpDy--list.r
# desc: Set scenario groups for Rmd rendering
# in: Rscript 20260825_cpDy--list.r
# out: upload/20260825_cpDy--list.csv
# arg: 0
# date: 20260826

b0 = list.files("../upload", full.names = T) 
hOst = read.csv(b0[grep("host", b0)], header = T)
sCe = read.csv(b0[grep("scenario", b0)], header = T)

a = hOst[unique(sCe$host),-(2:3)]

write.table(a, "../upload/20260825_cpDy--list.csv", row.names = F, col.names = F, sep = ",", quote = F)
