#!/bin/env Rscript
# author: ph-u
# script: analysis.r
# desc: self-defined functions for analyses on the transposon ecology simulations
# in: Rscript analysis.r
# out: NA
# arg: NA
# date: 20260703

source("func_analysis.r")
f0 = list.files("../data", ".rda", full.names = T)
pAr = read.csv("../raw/scenario.csv", header = T)

##### Load all data files #####
d0.host = d0.tpns = d0.offs = vector(mode = "list", length = length(f0))
names(d0.host) = names(d0.offs) = names(d0.tpns) = sub(".rda","",list.files("../data", ".rda", full.names = F))
for(i in 1:length(f0)){
  load(f0[i])
  d0.host[[i]] = rec.host
  d0.offs[[i]] = fAmily(rec.offspring)
  d0.tpns[[i]] = rec.transposon
};rm(i, rec.host, rec.offspring, rec.transposon)

##### overall data: Transposon dominance, Family tree #####
# https://biostatsquid.com/alpha-diversity-metrics/
rAtes = read.table(text = gsub("[-]","_",names(d0.tpns)), sep = "_")[,3:4]
colnames(rAtes) = c("seed","pAr")
rAtes$name = names(d0.tpns)
tpn.col = c("jump.rate","recombine.rate","generation","replicate","name","titre","shannon","simpson.reci","chao1", paste0("g", fam.typ <- unique(d0.offs[[1]][,1])[order(unique(d0.offs[[1]][,1]))]))
tpn.div = as.data.frame(matrix(0, nrow = length(d0.tpns)*nrow(d0.tpns[[1]]),  ncol = length(tpn.col)))
colnames(tpn.div) = tpn.col
tpn.div[,1] = rep(pAr$transposon[rAtes$pAr], each = nrow(d0.tpns[[1]]))
tpn.div[,2] = rep(pAr$gene[rAtes$pAr], each = nrow(d0.tpns[[1]]))
tpn.div[,3] = (1:nrow(d0.tpns[[1]]))-1
tpn.div[,4] = rep(rAtes$seed, each = nrow(d0.tpns[[1]]))
tpn.div[,5] = rep(rAtes$name, each = nrow(d0.tpns[[1]]))
for(i in 1:nrow(tpn.div)){
  tpn.i = lengths(strsplit(as.character(d0.tpns[[which(names(d0.tpns) == tpn.div$name[i])]][tpn.div$generation[i]+1,]), ";"))
  fam.i = table(as.numeric(d0.offs[[which(names(d0.tpns) == tpn.div$name[i])]][tpn.div$generation[i]+1,]))
  tpn.div[i,6:ncol(tpn.div)] = c(sum(tpn.i), diversity(tpn.i, index = "shannon"), diversity(tpn.i, index = "inv"), chao1(tpn.i, taxa.row = F),fam.i[match(fam.typ, names(fam.i))]/sum(fam.i))
};rm(i, tpn.i, fam.i)
tpn.div[is.na(tpn.div)] = 0
write.csv(tpn.div, "../data/tpn_dyn.csv", row.names = F, quote = F)
