#!/bin/env Rscript
# author: ph-u
# script: setScenario.r
# desc: set scenario into a csv file
# in: Rscript setScenario.r
# out: raw/scenario.csv
# arg: 0
# date: 20260707,20260721

##### Set simulation scenarios #####
tPn = c(10^-(1:4),0)
jpH1 = c("fixed", "charlesworth", "evolving")
pAr = list(
  copyRate = c(10^-(1:4),0),
  copy = c("fixed", "charlesworth", "evolving"),
  copyDir = c("both", "terminus", "origin"),
  gene = c(10^-(1:4),0),
  recom = c("switch", "homeostatic"),
  genotoxic = 0:4 # number of genotoxic events
)

a = data.frame(jumpRate = rep(tPn, each = length(jpH1)), jump = rep(jpH1, length(tPn)))
for(i in 1:length(pAr)){
  a = cbind(a, rep(pAr[[i]], each = nrow(a)))
};rm(i)

colnames(a)[-(1:2)] = names(pAr)
write.csv(a, "../raw/scenario.csv", row.names = F, quote = F)
