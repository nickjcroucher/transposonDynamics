#!/bin/env Rscript
# author: ph-u
# script: sum_sort.r
# desc: sort summary statistics
# in: bsub < sum_sort.r
# out: data/sum_sort.rda
# arg: 0
# date: 20260825

#BSUB -G team377f
#BSUB -o ../work/pj02-%J-%I.o
#BSUB -e ../work/pj02-%J-%I.e
#BSUB -q normal
#BSUB -M 20000
#BSUB -R "select[mem>20000] rusage[mem=20000] span[hosts=1]"
#BSUB -J "sum"

##### env #####
cat(date(),": set assemble env\n")
f0 = list.files("../data", pattern = paste0("sum_.*\\.rda$"), full.names = T)
b0 = list.files("../upload", full.names = T)
f.in = read.csv(b0[grep("input", b0)], header = T)
maxGen = as.numeric(f.in$Value[f.in$Type=="host organism constant generation number"])
tpn = read.csv(b0[grep("tpn", b0)], header = T)
hOst = read.csv(b0[grep("host", b0)], header = T)
sCe = read.csv(b0[grep("scenario", b0)], header = T)

gEt0 = cbind(tpn[match(sCe$transposon, tpn$uniqID),-(1:4)], hOst[match(sCe$host, row.names(hOst)),], f0[match(row.names(sCe), read.table(text = gsub("[.]rd", "@", gsub("_", "@", f0)), sep = "@")[,2])])
colnames(gEt0)[ncol(gEt0)] = "fNam"
row.names(gEt0) = NULL
gEt0 = cbind(gEt0, gen = rep(seq_len(maxGen+1)-1, each = nrow(gEt0)))
gEt0 = gEt0[order(gEt0$fNam),]

f0 = unique(gEt0$fNam)

##### Assemble plot df #####
cat(date(),": assemble plot df\n")
for(i in seq_len(length(f0))){
  cat(date(),":", i, "/", length(f0), "(",round(i/length(f0)*100,2),"% )\n")
  load(f0[i])
  t.hGene1 = t(t.hGene1)
  t.hGene2 = t(t.hGene2)
  row.names(t.hGene1) = row.names(t.hGene2) = NULL
  if(i > 1){
    fT0 = rbind(fT0, famTree)
    hSw = rbind(hSw, h.sweep)
    tGen = rbind(tGen, t.gen)
    tGir = rbind(tGir, t.gir)
    tPp = rbind(tPp, t.pop)
    tAv = rbind(tAv, t.ppavg)
    tH1 = merge(tH1, t.hGene1, all = T)
    tH2 = merge(tH2, t.hGene2, all = T)
  }else{
    fT0 = famTree
    hSw = h.sweep
    tGen = t.gen
    tGir = t.gir
    tPp = t.pop
    tAv = t.ppavg
    tH1 = t.hGene1
    tH2 = t.hGene2
  }
};rm(i)

##### export #####
save(gEt0, fT0, hSw, tGen, tGir, tPp, tAv, tH1, tH2, file = "../data/sum_sort.rda", compress = "xz")
