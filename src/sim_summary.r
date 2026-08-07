#!/bin/env Rscript
# author: ph-u
# script: sim_summary.r
# desc: simulation prelimiary summary
# in: source("sim_summary.r")
# out: NA
# date: 20260807

##### extra env #####
for(i in c("ape", "stringdist")){library(i, character.only=T)};rm(i)

##### Run analysis #####
cat(date(),": analysis started",argv[3],"-",argv[5],"\n")
gEnealogy = fAmily(rec.offspring)/ncol(rec.offspring)
rec.tpn = tpn.gen(rec.transposon)

##### Host organism domination #####
cat(date(),": map initial host organism dominance\n")
r.os0 = rec.offspring
r.os0[1,] = paste0(1:ncol(r.os0),";",1:ncol(r.os0))
rec.hostTraject = fAmily(r.os0)/ncol(r.os0)

##### Transposon perspective map #####
cat(date(),": map transposon dominance\n")
a.tpn = unique(unlist(strsplit(unlist(rec.transposon), ";")))
rec.tpnTraject = cbind(a.tpn, as.data.frame(matrix(0, nrow = length(a.tpn), ncol = nrow(rec.transposon))))
colnames(rec.tpnTraject) = c("transposon", paste0("gen",(1:nrow(rec.transposon))-1))
for(i in 1:nrow(rec.transposon)){
  i0 = table(unlist(strsplit(as.character(rec.transposon[i,]), ";")))
  for(i1 in 1:length(i0)){
    rec.tpnTraject[which(rec.tpnTraject[,1]==names(i0)[i1]),i+1] = i0[i1]
  };rm(i1)
};rm(i, i0)

##### Gene insertion map #####
cat(date(),": map recipient genes\n")
a.tpn.df = cbind(tPn.io(paste0(a.tpn, collapse = ";")), a.tpn)
rec.tpnGene.df = unique(substr(a.tpn.df$gene,3,nchar(a.tpn.df$gene)))
rec.tpnGene.df = cbind(rec.tpnGene.df, as.data.frame(matrix(0, nrow = length(rec.tpnGene.df), ncol = nrow(rec.transposon))))
colnames(rec.tpnGene.df) = c("gene", paste0("gen",(1:nrow(rec.transposon))-1))
a.tpn.df$map = match(substr(a.tpn.df$gene,3,nchar(a.tpn.df$gene)), rec.tpnGene.df$gene)
for(i in 1:nrow(rec.tpnGene.df)){
  rec.tpnGene.df[i,-1] = colSums(rec.tpnTraject[which(a.tpn.df$map==i),-1])
};rm(i)

##### Transposon tag map #####
cat(date(),": map transposon type dominance\n")
if(length(table(a.tpn.df$uniqID))>1){
  rec.uniqID = unique(a.tpn.df$uniqID)
  rec.uniqID = cbind(rec.uniqID, as.data.frame(matrix(0, nrow = length(rec.uniqID), ncol = nrow(rec.transposon))))
  colnames(rec.uniqID) = c("uniqID", paste0("gen",(1:nrow(rec.transposon))-1))
  a.tpn.df$uMap = match(a.tpn.df$uniqID, rec.uniqID)
  for(i in 1:nrow(rec.uniqID)){
    rec.uniqID[i,-1] = colSums(rec.tpnTraject[which(a.tpn.df$uMap==i),-1])
  };rm(i)
}else{
  rec.uniqID = cbind(names(table(a.tpn.df$uniqID)), colSums(rec.tpnTraject[,-1]))
}

##### Host genome phylogenetics ##### !!!

##### Export #####
save(gEnealogy, rec.tpn, rec.hostTraject, rec.tpnTraject, rec.tpnGene.df, rec.uniqID, file = paste0("../data/ana--", argv[3], "_", argv[5], ".rda"), compress = "xz")
cat(date(),": simulation and analysis completed",argv[3],"-",argv[5],"\n")
