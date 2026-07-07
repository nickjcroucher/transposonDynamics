#!/bin/env Rscript
# author: ph-u
# script: setScenario.r
# desc: set scenario into a csv file
# in: Rscript setScenario.r
# out: raw/scenario.csv
# arg: 0
# date: 20260707

tPn = 10^-(0:5)
gEne = 10^-(0:5)
jpH1 = c("fixed", "charlesworth", "evolving", "genotoxic", "evolving_genotoxic")
grH1 = c("croucher", "adv_croucher")

a = data.frame(transposon = rep(tPn, each = length(gEne)), gene = rep(gEne, length(tPn)))
a0 = data.frame(transposon = rep(a[,1], length(jpH1)), gene = rep(a[,2], length(jpH1)), jump = rep(jpH1, each = nrow(a)))
a = data.frame(transposon = rep(a0[,1], length(grH1)), gene = rep(a0[,2], length(grH1)), jump = rep(a0[,3], length(grH1)), recom = rep(grH1, each = nrow(a0)))

write.csv(a, "../raw/scenario.csv", row.names = F, quote = F)
