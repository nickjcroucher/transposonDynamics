#!/bin/env Rscript
# author: ph-u
# script: 20260904_tpnEqm-comHist.r
# desc: Compare transposon spread on gene recombination mechanisms between survived and eliminated populations
# in: Rscript 20260904_tpnEqm-comHist.r
# out: upload/20260904_tpnEqm-comHist.pdf
# arg: 0
# date: 20260907

load("../data/20260904_tpnEqm.rda")
s.sum = g.Hit0
s.all = g.Hit
s.all[,1] = as.numeric(as.character(s.all[,1]))
s.tab = table(s.all[,1])
print(s.tab[length(s.tab)]/as.numeric(names(s.tab)[length(s.tab)]))

load("../data/20260904_tpnEqm-elim.rda")
d.sum = g.Hit0
d.all = g.Hit
d.all[,1] = as.numeric(as.character(d.all[,1]))
d.tab = table(d.all[,1])
print(d.tab[length(d.tab)]/as.numeric(names(d.tab)[length(d.tab)]))
rm(g.Hit,g.Hit0)

##### Transposon-hitting recombination mechanism #####
## higher chance to be eliminated ultimately if dropping at first dozens of generations
## higher chance to survive ultimately if increase at first dozens of generations (significant)
w0 = data.frame(genHit = unique(s.sum$genHit), W = NA, p = NA)
for(i in seq_len(nrow(w0))){
  i0 = wilcox.test(d.all$Freq[which(d.all$genHit==w0[i,1])], s.all$Freq[which(s.all$genHit==w0[i,1])])
  w0[i,-1] = c(i0$statistic, i0$p.value)
};rm(i,i0)
w0$p.adj = p.adjust(w0$p, method = "BH")
w0[w0$p.adj<.05,]

#pdf("../upload/20260904_tpnEqm-comHist.pdf")
pdf("../res/20260904_tpnEqm-comHist.pdf", width = 10, height = 7)
par(mar = c(5,4,4,0)+.1, mfrow = c(2,2))

hist(d.all$Freq[which(d.all$genHit==0)], col = "#00000033", probability=T, breaks=100, border = "#00000000", xlim = c(0,.1),
  main = "Grey = !0.1 gene recombination rate\nBlue = 0.1 (gen = 0)", xlab = "Proportion of population with transposons\nhitting gene recombination mechanism")
hist(s.all$Freq[which(s.all$genHit==0)], col = "#0000ff33", probability=T, breaks=100, border = "#00000000", add=T)

hist(d.all$Freq[which(d.all$genHit > 0 & d.all$genHit < 10)], col = "#00000033", probability=T, breaks=100, border = "#00000000", xlim = c(0,1),
  main = "Grey = !0.1 gene recombination rate\nBlue = 0.1 (0 < gen < 10)", xlab = "Proportion of population with transposons\nhitting gene recombination mechanism")
hist(s.all$Freq[which(s.all$genHit > 0 & s.all$genHit < 10)], col = "#0000ff33", probability=T, breaks=100, border = "#00000000", add=T)

hist(d.all$Freq[which(d.all$genHit > 10 & d.all$genHit < 500)], col = "#00000033", probability=T, breaks=100, border = "#00000000", xlim = c(0,1),
  main = "Grey = !0.1 gene recombination rate\nBlue = 0.1 (10 < gen < 500)", xlab = "Proportion of population with transposons\nhitting gene recombination mechanism")
hist(s.all$Freq[which(s.all$genHit > 10 & s.all$genHit < 500)], col = "#0000ff33", probability=T, breaks=100, border = "#00000000", add=T)

hist(d.all$Freq[which(d.all$genHit==500)], col = "#00000033", probability=T, breaks=100, border = "#00000000", xlim = c(0,1),
  main = "Grey = !0.1 gene recombination rate\nBlue = 0.1 (gen = 500)", xlab = "Proportion of population with transposons\nhitting gene recombination mechanism")
hist(s.all$Freq[which(s.all$genHit==500)], col = "#0000ff33", probability=T, breaks=100, border = "#00000000", add=T)

invisible(dev.off())
