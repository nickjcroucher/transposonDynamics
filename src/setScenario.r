#!/bin/env Rscript
# author: ph-u
# script: setScenario.r
# desc: set scenario into a csv file
# in: Rscript setScenario.r
# out: raw/scenario.csv
# arg: 0
# date: 20260707,20260721

set.seed(1234)

##### Set transposon simulation scenarios #####
gEne = ""
lOc = 0
pAr = list(
  generation = 0,
  valid = T,
  uniqID = "",
  size = c(500,1000,1500,2000),
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
uID.len = ceiling(log(nrow(a))/log(length(LETTERS)))
repeat{
  uID = unique(as.data.frame(matrix(sample(LETTERS, nrow(a)*(uID.len+1), replace = T), ncol = uID.len+1)))
  if(nrow(uID) >= nrow(a)){break}
}
a$uniqID = apply(uID,1,paste0, collapse = "")[1:nrow(a)]

write.csv(a, "../raw/template-tpn.csv", row.names = F, quote = F)

##### Set cell simulation scenarios #####
rEcom = c(10^-(1:4),0)
rEcH1 = c("switch", "homeostatic")
pAr = list(
  cell = c("haploid", "diploid"),
  transposonEffect = c(T,F),
  genotoxic = 0:4 # number of genotoxic events
)

a0 = data.frame(recom = rep(rEcom, each = length(rEcH1)), recomH1 = rep(rEcH1, length(rEcom)))
for(i in 1:length(pAr)){
  a0 = cbind(a0, rep(pAr[[i]], each = nrow(a0)))
};rm(i)

colnames(a0)[-(1:2)] = names(pAr)
write.csv(a0, "../raw/template-host.csv", row.names = F, quote = F)

##### Set overall scenario to do simulation #####
tPn = a$uniqID[which(a$size == 1000 & a$jumpH1 == "fixed" & a$copyH1 == "fixed" & a$copyDir == "both")]
hOst = row.names(a0)[which(a0$recomH1 == "switch" & a0$cell == "haploid")]

res = data.frame(transposon = rep(tPn, each = length(hOst)), host = hOst)
write.csv(res, "../raw/scenario.csv", row.names = F, quote = F)
