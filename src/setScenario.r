#!/bin/env Rscript
# author: ph-u
# script: setScenario.r
# desc: set scenario into a csv file
# in: Rscript setScenario.r
# out: raw/scenario.csv
# arg: 0
# date: 20260707,20260721

##### Set transposon simulation scenarios #####
cat(date(),": set transposon table\n")
gEne = ""
lOc = 0
pAr = list(
  generation = 0,
  valid = T,
  uniqID = "",
  size = 1000,
  jumpRate = c(10^-(1:4),0),
  jumpH1 = c("fixed", "charlesworth", "evolving"),
  copyRate = c(10^-(1:4),0),
  copyH1 = c("fixed", "charlesworth", "evolving"),
  copyDir = c("both", "terminus", "origin")
)

a = data.frame(gene = rep(gEne, each = length(lOc)), location = rep(lOc, length(gEne)))
for(i in 1:length(pAr)){
  a = cbind(a, rep(pAr[[i]], each = nrow(a)))
};rm(i)

colnames(a)[-(1:2)] = names(pAr)

## Set transposon uniqID
cat(date(),": set transposon tags\n")
uID.len = ceiling(log(nrow(a))/log(length(LETTERS)))
uID.set = as.data.frame(cbind(LETTERS, rep(LETTERS, each = length(LETTERS))))
if(uID.len > 2){for(i in 3:uID.len){
   uID.set = cbind(uID.set, rep(LETTERS, each = nrow(uID.set)))
};rm(i)}
a$uniqID = apply(uID.set,1,paste0, collapse = "")[1:nrow(a)]

write.csv(a, "../raw/template-tpn.csv", row.names = F, quote = F)

##### Set cell simulation scenarios #####
cat(date(),": set host cell table\n")
rEcom = c(10^-(1:4),0)
rEcH1 = c("switch", "homeostatic")
pAr = list(
  cell = c("haploid", "diploid"),
  homologousAutoRecom = (seq_len(11)-1)/10,
  transposonEffect = c(T,F),
  genotoxic = 0:4 # number of genotoxic events
)

a0 = data.frame(recom = rep(rEcom, each = length(rEcH1)), recomH1 = rep(rEcH1, length(rEcom)))
for(i in 1:length(pAr)){
  a0 = cbind(a0, rep(pAr[[i]], each = nrow(a0)))
};rm(i)

colnames(a0)[-(1:2)] = names(pAr)
row.names(a0) = NULL
write.csv(a0, "../raw/template-host.csv", row.names = F, quote = F)

##### Set overall scenario to do simulation #####
cat(date(),": set scenario table\n")
sCe = data.frame(transposon = a$uniqID, host = rep(row.names(a0), each = nrow(a)))
write.csv(sCe, "../raw/scenario.csv", row.names = F, quote = F)

##### Set rerun scenarios #####
tPn = a$uniqID[which(a$jumpH1 == "fixed" & a$copyH1 == "fixed" & a$copyDir == "both")]
hOst = row.names(a0)[which(a0$recom == 0 & a0$recomH1 == "switch" & a0$cell == "haploid" & a0$homologousAutoRecom %in% c(0, .1, .5, .9, 1) & a0$transposonEffect == F & a0$genotoxic == 0)]

res = which(sCe$transposon %in% tPn & sCe$host %in% hOst) #data.frame(transposon = rep(tPn, each = length(hOst)), host = hOst)
res = as.data.frame(cbind(rep(seq_len(10), each = length(res)), res))
write.table(res, "../data/rerun.csv", sep = ",", row.names = F, col.names = F, quote = F)
