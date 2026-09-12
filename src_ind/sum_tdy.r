#!/bin/env Rscript
# author: ph-u
# script: sum_tdy.r
# desc: summarize pj02 raw simulation output
# in: Rscript sum_tdy.r
# out: sum_tdy--[simulation set].rda
# arg: 1
# date: 20260824

argv = as.numeric(commandArgs(T))
if(length(argv)==0){argv = 1}
f.sce = read.csv("../upload/scenario.csv", header = T)
f.rda = list.files("../data", pattern = paste0("ana--.*_",argv,"\\.rda$"), full.names = T)

##### Calculate averages #####
for(i in seq_len(length(f.rda))){
  load(f.rda[i])

  rec.hostTraject = rec.hostTraject[,rev(order(as.numeric(rec.hostTraject[nrow(rec.hostTraject),])))]
  colnames(rec.hostTraject) = paste0("h",seq_len(ncol(rec.hostTraject)))

  if(i > 1){
    famTree = famTree + gEnealogy
    h.sweep = h.sweep + rec.hostTraject
    t.pop = t.pop + as.numeric(rec.uniqID[,2])
    t.gir = t.gir + rec.tpn$distribution$g/rowSums(rec.tpn$distribution)
    t.ppavg = t.ppavg + rowMeans(rec.tpn$count)
    t.gen = t.gen + rec.tpn$generation

    t.hGene = rbind(t.hGene, rec.tpnGene.df)
  }else{
    famTree = gEnealogy
    h.sweep = rec.hostTraject
    t.pop = as.numeric(rec.uniqID[,2])
    t.gir = rec.tpn$distribution$g/rowSums(rec.tpn$distribution)
    t.ppavg = rowMeans(rec.tpn$count)
    t.gen = rec.tpn$generation

    t.hGene = rec.tpnGene.df
  }
};rm(i)

famTree = famTree / length(f.rda)
h.sweep = h.sweep / length(f.rda)
t.pop = t.pop / length(f.rda)
t.gir = t.gir / length(f.rda)
t.ppavg = t.ppavg / length(f.rda)
t.gen = t.gen / length(f.rda)
t.hGene1 = (t.hGene2 <- rowsum(t.hGene[,-1], t.hGene[,1]))/as.vector(table(t.hGene[,1]))

save(famTree, h.sweep, t.pop, t.gir, t.ppavg, t.gen, t.hGene1, t.hGene2, file = paste0("../data/sum_",argv,".rda"))

