#!/bin/env Rscript
# author: ph-u
# script: ai_convergence.r
# desc: paths
# in: source("ai_convergence.r")
# out: NA
# arg: 0
# date: 20261005

source("src.r")
library(jsonlite)

table(f.0[,-1])
f.ref[f.ref$Accession %in% f.0$Accession[f.0[,2]=="Transposase" & f.0[,3]=="Unknown"],]
f.ref[f.ref$Accession %in% f.0$Accession[f.0[,2]=="Unknown" & f.0[,3]=="Cellular functions"],]

f.1 = f.0[which(f.0[,2] %in% c("Reverse transcriptase", "Transposase") | f.0[,3] %in% c("Reverse transcriptase", "Transposase")),]
write.table(f.1, "../data/CL0219_category.csv", sep = ",", col.names = F, quote = F, row.names = F)

##### Env: CL0219 #####
d.0 = as.data.frame(fromJSON("../raw/CL0219_raw.jsonl"))
write.table(d.0$results.metadata$accession, "../raw/CL0219_tax_report_raw.txt", quote = F, row.names = F, col.names = F)
#d.0 = cbind(d.0$results.metadata, d.0$results.sets) # infinitely loading
colnames(d.0$results.metadata)
d.1 = as.data.frame(t(matrix(unlist(d.0$results.sets), nrow = 3)))
system("sed -e 's/|\t/@/g' -e 's/\t//g' -e 's/ @/@/g' ../raw/CL0219_tax_report.txt > ../raw/CL0219_tax_report-1.txt") # excel delimiter text to columns
d.2 = read.table("../raw/CL0219_tax_report-1.txt", header = T, sep = "\t", quote = "", comment.char = "")
write.table(unique(c(d.2$taxid, unlist(strsplit(d.2$lineage, " ")))), "../raw/CL0219_tax_report_rawFull.txt", quote = F, row.names = F, col.names = F)
system("sed -e 's/|\t/@/g' -e 's/\t//g' -e 's/ @/@/g' ../raw/CL0219_tax_report-Full.txt > ../raw/CL0219_tax_report-Full-1.txt") # excel delimiter text to columns
d.2 = read.table("../raw/CL0219_tax_report-Full-1.txt", header = T, sep = "\t", quote = "", comment.char = "")

