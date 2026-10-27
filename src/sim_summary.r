#!/bin/env Rscript
# author: ph-u
# script: sim_summary.r
# desc: simulation prelimiary summary
# in: source("sim_summary.r")
# out: NA
# date: 20260807

##### extra env #####
source("func_analysis.r")
rec.gRecom0 = rec.gRecom[seq_len(gEn), , drop = FALSE]
rec.host0 = rec.host[seq_len(gEn), , drop = FALSE]
rec.offspring0 = rec.offspring[seq_len(gEn), , drop = FALSE]
rec.tpnCount0 = rec.tpnCount[seq_len(gEn), , drop = FALSE]
rec.transposon0 = rec.transposon[seq_len(gEn), , drop = FALSE]
# for(i in c("ape", "stringdist")){library(i, character.only=T)};rm(i)
pdf(paste0("../data/ana--",argv[3],"-",argv[5],".pdf"))

##### Host parameters mapping #####
# cat(date(),": Genotype dominance",argv[3],"-",argv[5],"\n")
gEnealogy = fAmily(rec.offspring0, host.0$cell)/ncol(rec.offspring0)
suppressWarnings(matplot(x = seq_len(nrow(gEnealogy))-1, y = gEnealogy*100, col = cBp, xlim = c(0, gEn.max), ylim = c(0,100), xlab = "Generation", ylab = "Dominance (%)", main = paste0("Genealogy of genotypes (df = gEnealogy)\nrep = ",argv[3],", scenario = ",argv[5])))

# cat(date(),": Initial host dominance",argv[3],"-",argv[5],"\n")
r.os0 = rec.offspring0
r.os0[1,] = paste0(seq_len(ncol(r.os0)),";",seq_len(ncol(r.os0)))
rec.hostTraject = fAmily(r.os0, "haploid")/ncol(r.os0)
suppressWarnings(matplot(x = seq_len(nrow(rec.hostTraject))-1, y = rec.hostTraject*100, col = cBp, xlim = c(0, gEn.max), ylim = c(0,100), xlab = "Generation", ylab = "Dominance (%)", main = paste0("Genealogy of Initial host (df = rec.hostTraject)\nrep = ",argv[3],", scenario = ",argv[5])))

##### Transposon parameters mapping #####
# cat(date(),": Transposon counts",argv[3],"-",argv[5],"\n")
plot(x = seq_len(gEn)-1, y = log10(rowSums(rec.tpnCount0)+1), xlim = c(0, gEn.max), ylim = c(0, max(log10(rowSums(rec.tpnCount0)+1))), xlab = "Generation", ylab = "Total transposon copy number [log10(y+1)]", main = paste0("Total transposon copy number [N = ",hPop,"] \n(df = rec.tpnCount) rep = ",argv[3],", scenario = ",argv[5]))
plot(x = seq_len(gEn)-1, y = rowSums(rec.tpnCount0) / hPop, xlim = c(0, gEn.max), ylim = c(0, max(rowSums(rec.tpnCount0)) / hPop), xlab = "Generation", ylab = "Average transposon copy number", main = paste0("Average transposon copy number [N = ",hPop,"] \n(df = rec.tpnCount) rep = ",argv[3],", scenario = ",argv[5]))

# cat(date(),": Transposon activity dominance",argv[3],"-",argv[5],"\n")
tpn.dom = matrix(0, nrow = gEn, ncol = gEn)

gNam = unique(substr(geneIn$gene$locus_tag,2,nchar(geneIn$gene$locus_tag)))
tpn.gene = matrix(0, nrow = gEn, ncol = length(gNam))
colnames(tpn.gene) = gNam

g.Fit = matrix(0, nrow = gEn, ncol = 10)
g.Rec = matrix(0, nrow = gEn, ncol = 2)
g.Fit.ref = geneIn$gene[which(substr(geneIn$gene$locus_tag,1,1)==0),c("locus_tag", "fitness", "recombination")]
g.Fit.ref$locus_tag = substr(g.Fit.ref$locus_tag,2,nchar(g.Fit.ref$locus_tag))
g.Fit.ref$fitness = round(g.Fit.ref$fitness,1)
g.Fit.ref$recombination = ifelse(is.na(g.Fit.ref$recombination),0,1)

for(i in seq_len(nrow(rec.transposon0))){
  t.tmp = tPn.io(as.character(rec.transposon0[i,]))

  ## Generation dominance
  t.gen = table(t.tmp$generation)
  t.gen = t.gen/sum(t.gen)
  for(i0 in seq_len(length(t.gen))){tpn.dom[i,as.numeric(names(t.gen))+1] = t.gen};rm(i0)

  ## Gene dominance
  t.gen = table(substr(t.tmp$gene,3,nchar(t.tmp$gene)))
  for(i0 in seq_len(length(t.gen))){tpn.gene[i,names(t.gen)] = t.gen};rm(i0)
};rm(i, t.tmp, t.gen)

## Gene essentiality
for(i in seq_len(ncol(g.Fit))){
  g.Fit[,i] = rowSums(tpn.gene[,colnames(tpn.gene) %in% g.Fit.ref$locus_tag[g.Fit.ref$fitness==i/10], drop = FALSE])
};rm(i)
g.Fit = g.Fit/rowSums(g.Fit)

## Gene recombination mechanism
for(i in seq_len(ncol(g.Rec))){
  g.Rec[,i] = rowSums(tpn.gene[,colnames(tpn.gene) %in% g.Fit.ref$locus_tag[g.Fit.ref$recombination==(i-1)], drop = FALSE])
};rm(i)
g.Rec[,1] = log1p(g.Rec[,1]); g.Rec = g.Rec[,c(2,1), drop = FALSE]

suppressWarnings(matplot(x = seq_len(nrow(tpn.dom))-1, y = tpn.dom*100, col = cBp, xlim = c(0, gEn.max), ylim = c(0,100), xlab = "Generation", ylab = "Dominance (%)", main = paste0("Age of transposons' last known activity (df = tpn.dom)\nrep = ",argv[3],", scenario = ",argv[5])))
suppressWarnings(matplot(x = seq_len(nrow(g.Fit))-1, y = g.Fit*100, col = cBp, xlim = c(0, gEn.max), ylim = c(0,100), xlab = "Generation", ylab = "Dominance (%)", main = paste0("Fitness of genes hosting transposons (df = g.Fit from tpn.gene)\nrep = ",argv[3],", scenario = ",argv[5])))
suppressWarnings(matplot(x = seq_len(nrow(g.Rec))-1, y = g.Rec, col = cBp, xlim = c(0, gEn.max), xlab = "Generation", ylab = "Count", main = paste0("Count of transposons hitting gene recombination \nmechanism and those which are not [log1p] (df = g.Rec)\nrep = ",argv[3],", scenario = ",argv[5])))

##### Export #####
invisible(dev.off())
save(gEnealogy, rec.hostTraject, rec.tpnCount, tpn.dom, tpn.gene, g.Fit, g.Rec, g.Fit.ref, file = paste0("../data/ana--", argv[3], "_", argv[5], ".rda"), compress = "xz")
