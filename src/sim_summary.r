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
# cat(date(),": map transposon dominance\n")
a.tpn = unique(unlist(strsplit(unlist(rec.transposon), ";")))
# rec.tpnTraject = cbind(a.tpn, as.data.frame(matrix(0, nrow = length(a.tpn), ncol = nrow(rec.transposon))))
# colnames(rec.tpnTraject) = c("transposon", paste0("gen",seq_len(nrow(rec.transposon))-1))
# tpnSplit = lapply(seq_len(nrow(rec.transposon)), function(i){
#   v = unlist(strsplit(unlist(rec.transposon[i,], use.names = F), ";"))
#   return(v[nzchar(v)])
# })
# a.tpn = unique(unlist(tpnSplit))
# tpnMat = matrix(0, nrow = length(a.tpn), ncol = length(tpnSplit))
# for(i in seq_along(tpnSplit)){
#   tpnMat[,i] = tabulate(match(tpnSplit[[i]], a.tpn), nbins = length(a.tpn))
# }
# rec.tpnTraject = data.frame(transposon = a.tpn, tpnMat, stringsAsFactors = F)
# colnames(rec.tpnTraject) = c("transposon", paste0("gen", seq_len(nrow(rec.transposon))-1))

##### Gene insertion map #####
cat(date(),": map recipient genes\n")
a.tpn.df = cbind(tPn.io(paste0(a.tpn, collapse = ";")), a.tpn)
rec.tpnGene.df = unique(substr(a.tpn.df$gene,3,nchar(a.tpn.df$gene)))
rec.tpnGene.df = cbind(rec.tpnGene.df, as.data.frame(matrix(0, nrow = length(rec.tpnGene.df), ncol = nrow(rec.transposon))))
colnames(rec.tpnGene.df) = c("gene", paste0("gen",(1:nrow(rec.transposon))-1))
a.tpn.df$map = match(substr(a.tpn.df$gene,3,nchar(a.tpn.df$gene)), rec.tpnGene.df$gene)
# tmp = rowsum(as.matrix(rec.tpnTraject[,-1]), group = a.tpn.df$map, reorder = T)
# idx = match(seq_len(nrow(rec.tpnGene.df)), as.integer(rownames(tmp)))
# rec.tpnGene.df[!is.na(idx), -1] = tmp[idx[!is.na(idx)], , drop = F]

##### Transposon tag map #####
# cat(date(),": map transposon type dominance\n")
# if(length(table(a.tpn.df$uniqID))>1){
#   rec.uniqID = unique(a.tpn.df$uniqID)
#   rec.uniqID = cbind(rec.uniqID, as.data.frame(matrix(0, nrow = length(rec.uniqID), ncol = nrow(rec.transposon))))
#   colnames(rec.uniqID) = c("uniqID", paste0("gen",(1:nrow(rec.transposon))-1))
#   a.tpn.df$uMap = match(a.tpn.df$uniqID, rec.uniqID)
#   for(i in 1:nrow(rec.uniqID)){
#     rec.uniqID[i,-1] = colSums(rec.tpnTraject[which(a.tpn.df$uMap==i),-1])
#   };rm(i)
# }else{
#   rec.uniqID = cbind(names(table(a.tpn.df$uniqID)), colSums(rec.tpnTraject[,-1]))
# }

##### Gene recombination mechanism hit ratio #####
hPos = grep(gsub(";","|",inFile$params$Value[inFile$params$Type=="genes for recombination mechanism"]), unlist(rec.transposon))
if(length(hPos) > 0){
  h.Gen = (hPos %% nrow(rec.transposon))-1
  h.Gen[h.Gen < 0] = nrow(rec.transposon)-1
  h.Gen = as.data.frame(table(h.Gen)/ncol(rec.transposon))
  h.Gen[,1] = as.numeric(as.character(h.Gen[,1]))
}else{
  h.Gen = data.frame(Var1 = seq_len(nrow(rec.transposon))-1, Freq = 0)
}
h.Miss = (seq_len(nrow(rec.transposon))-1)[!((seq_len(nrow(rec.transposon))-1) %in% h.Gen[,1])]
if(length(h.Miss) > 0){
  h.Miss = as.data.frame(matrix(c(h.Miss, rep(0, length(h.Miss))), nrow = length(h.Miss)))
  colnames(h.Miss) = colnames(h.Gen)
  h.Gen = rbind(h.Gen, h.Miss)
}
gRecom.hitRatio = h.Gen[order(h.Gen[,1]),]
rownames(gRecom.hitRatio) = NULL

##### Host genome phylogenetics ##### !!!

##### Export #####
save(gEnealogy, rec.tpn, rec.hostTraject, rec.tpnGene.df, gRecom.hitRatio, file = paste0("../data/ana--", argv[3], "_", argv[5], ".rda"), compress = "xz") #, rec.tpnTraject, rec.uniqID
