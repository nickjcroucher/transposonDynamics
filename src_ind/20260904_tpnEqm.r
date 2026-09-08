#!/bin/env Rscript
# author: ph-u
# script: 20260904_tpnEqm.r
# desc: proportion of recombination genes hit per generation in population
# in: bsub < 20260904_tpnEqm.r
# out: data/20260904_tpnEqm.rda
# arg: 0
# date: 20260904

# BSUB -G team377f
# BSUB -o ../work/pj02-%J-%I.o
# BSUB -e ../work/pj02-%J-%I.e
# BSUB -q normal
# BSUB -M 4000
# BSUB -R "select[mem>4000] rusage[mem=4000] span[hosts=1]"
# BSUB -J "tEq"

##### env #####
pLt = 1 # 0 or 1
source("colour.r")
f0 = list.files("../data", full.names = T)
b0 = list.files("../upload", full.names = T)
f.in = read.csv(b0[grep("input", b0)], header = T)
tpn = read.csv(b0[grep("tpn.csv", b0)], header = T)
hOst = read.csv(b0[grep("host", b0)], header = T)
sCe = read.csv(f0[grep("scenario", f0)], header = T)
rEr = read.csv(f0[grep("rerun", f0)], header = F)
gRecom = strsplit(f.in$Value[f.in$Type=="genes for recombination mechanism"], ";")[[1]]
maxGen = as.numeric(f.in$Value[f.in$Type=="host organism constant generation number"])
popNum = as.numeric(f.in$Value[f.in$Type=="host organism constant population size"])

hOst0 = hOst[which(hOst$recom==.1 & hOst$recomH1=="switch" & hOst$cell=="haploid"),]
tpn0 = tpn[which(tpn$jumpH1=="fixed" & tpn$copyH1=="fixed" & tpn$size==1000 & tpn$copyDir=="both"),]
sCe0 = sCe[which(sCe$host==row.names(hOst0) & sCe$transposon %in% tpn0$uniqID),]
if(pLt > 0){
  rUns = rEr[which(rEr$V2 %in% row.names(sCe0)),]
}else{
  rUns = rEr[sample(which(!(rEr$V2 %in% row.names(sCe0))), 500, replace = F),]
}
rEs = as.data.frame(matrix(0, nrow = maxGen, ncol = nrow(rUns)))

##### Recombination mechanism hit ratio #####
cat(date(),": recom mechanism hit ratio calculation\n")
for(i in seq_len(nrow(rUns))){
  cat(date(),":",i,"/",nrow(rUns),"(",round(i/nrow(rUns)*100,2),"% )\n")
  load(f0[grep(paste0("tPn--",paste0(rUns[i,], collapse = "_"), ".rda"), f0)])
  genHit = (grep(paste0(gRecom, collapse="|"), unlist(rec.transposon)) %% (maxGen+1))-1 # generation
  genHit[genHit < 0] = maxGen
  #ceiling(gHit / (maxGen+1)) # organism
  if(i>1){
      g.Hit = rbind(g.Hit, as.data.frame(table(genHit)/popNum))
  }else{
      g.Hit = as.data.frame(table(genHit)/popNum)
  }
};rm(i);cat(date(),": done\n")

##### Summary statistics #####
g.Hit0 = aggregate(Freq ~ genHit, data = g.Hit, FUN = quantile, probs = c(.05, .5, .95), simplify = T)

if(pLt > 0){
  save(g.Hit, g.Hit0, file = "../data/20260904_tpnEqm.rda", compress = "xz")
  pdf("../upload/20260904_tpnEqm.pdf", width = 21, height = 10)
}else{
  save(g.Hit, g.Hit0, file = "../data/20260904_tpnEqm-elim.rda", compress = "xz")
  pdf("../upload/20260904_tpnEqm-elim.pdf", width = 21, height = 10)
}
par(mar = c(5,5,0,0)+.1, cex = 2)
plot(x = g.Hit0$genHit, y = g.Hit0$Freq[,2], ylim = c(0, max(g.Hit0$Freq)), xlab = "Generation", ylab = "Individuals with gene recombination\nmechanism hit by transposons (%)")
polygon(x = c(as.numeric(as.character(g.Hit0$genHit)), rev(as.numeric(as.character(g.Hit0$genHit)))), y = c(g.Hit0$Freq[,1], rev(g.Hit0$Freq[,3])), col = "#00000040", border = NA)
invisible(dev.off())
