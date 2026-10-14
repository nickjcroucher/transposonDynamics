#!/bin/env Rscript
# author: ph-u
# script: ai_convergence.r
# desc: paths
# in: source("ai_convergence.r")
# out: NA
# arg: 0
# date: 20261005

source("src.r")

table(f.0[,-1])
f.ref[f.ref$Accession %in% f.0$Accession[f.0[,2]=="Transposase" & f.0[,3]=="Unknown"],]
f.ref[f.ref$Accession %in% f.0$Accession[f.0[,2]=="Unknown" & f.0[,3]=="Cellular functions"],]
