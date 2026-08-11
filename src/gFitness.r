#!/bin/env Rscript
# author: ph-u
# script: gFitness.r
# desc: quantify gene fitness
# in: Rscript gFitness.r [../path/2/input.csv] [../path/2/tnseq.csv]
# out: data/gFitness.r
# arg: 1
# date: 20260811

argv = (commandArgs(T))
if(length(argv) < 2){argv=c("../raw/input.csv", "../data/tnseq.csv")}

source("func.r")
library(glmmTMB)
inFile = inParams(argv[1])

a = read.csv(argv[2], header = T)
a$length = inFile$gene$length[match(a$locus_tag, substr(inFile$gene$locus_tag, 2, nchar(inFile$gene$locus_tag)))]
a$density = a$count/a$length
g0 = glmmTMB(count ~ length + locus_tag, data = a, family = poisson(link = "identity"), ziformula = ~1)
g0 = as.data.frame(summary(g0)$coefficient$cond)
g0$gene = gsub("locus_tag","",row.names(g0))
g0$gene[1] = a$locus_tag[1]
a$fitness = g0$Estimate[match(a$locus_tag, g0$gene)]
a$fitness = a$fitness - a$fitness[1]
a$fitness = reZero(a$fitness, new0 = min(a$fitness, na.rm = T), new1 = max(a$fitness, na.rm = T))*2
a$fitness[is.na(a$fitness)] = 1
a$fitness[a$count==0] = 0

write.csv(a[,c("locus_tag", "fitness")], "../data/gFitness.csv", row.names = F, quote = F)
