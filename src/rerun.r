#!/bin/env Rscript
# author: ph-u
# script: rerun.r
# desc: set rerun scenarios and replicates
# in: Rscript rerun.r
# out: data/rerun.csv
# arg: 0
# date: 20260710

source("p_path.r")
f0 = read.csv(p.file[grep("scenario.csv", p.file)[1]], header = T)
#f0 = read.table(p.file[grep("yet", p.file)], header = F, sep = "_")
f0 = f0[which(f0$jump=="fixed" & f0$copy=="fixed" & f0$recom=="switch" & f0$copyDir=="both" & f0$genotoxic==1 & (f0$jumpRate==0 | f0$copyRate==0 | f0$gene==0)),]
d0 = data.frame(rep = rep(1:10, each = nrow(f0)), sce = 1:nrow(f0))
write.table(d0, "../data/rerun.csv", row.names = F, col.names = F, sep = ",")
