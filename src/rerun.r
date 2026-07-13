#!/bin/env Rscript
# author: ph-u
# script: rerun.r
# desc: set rerun scenarios and replicates
# in: Rscript rerun.r
# out: data/rerun.csv
# arg: 0
# date: 20260710

source("p_path.r")
f0 = read.table(p.file[grep("yet", p.file)], header = F, sep = "_")
write.table(f0, "../data/rerun.csv", row.names = F, col.names = F, sep = ",")
