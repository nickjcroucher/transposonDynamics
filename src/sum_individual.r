#!/bin/env Rscript
# author: ph-u
# script sum_individual.r
# desc: summary from the perspective of individuals
# in: Rscript sum_individual.r
# out: data/sum_ind.csv
# arg: 0
# date: 20260719

source("p_path.r")
source("func_analysis.r")
f.in = p.file[grep("ana_ind", p.file)]
f.in = f.in[grep("rda$", f.in)]
f.in = f.in[grep("Trans", f.in)]

##### Concatenate data files #####
for(i in 1:length(f.in)){
  load(f.in[i])
  r.gAll[is.na(r.gAll)] = r.jAll[is.na(r.jAll)] = r.nAll[is.na(r.nAll)] = 0
  for(i0 in 1:length(r.gGen)){
    r.gGen[[i0]][is.na(r.gGen[[i0]])] = r.jGen[[i0]][is.na(r.jGen[[i0]])] = r.nGen[[i0]][is.na(r.nGen[[i0]])] = 0
  }
  if(i>1){
    r0.gAll = r.gAll + r0.gAll
    r0.jAll = r.jAll + r0.jAll
    r0.nAll = r.nAll + r0.nAll
    for(i0 in 1:length(r.gGen)){
      r0.gGen[[i0]] = r.gGen[[i0]] + r0.gGen[[i0]]
      r0.jGen[[i0]] = r.jGen[[i0]] + r0.jGen[[i0]]
      r0.nGen[[i0]] = r.nGen[[i0]] + r0.nGen[[i0]]
    }
  }else{
    r0.gAll = r.gAll
    r0.jAll = r.jAll
    r0.nAll = r.nAll
    r0.gGen = r.gGen
    r0.jGen = r.jGen
    r0.nGen = r.nGen
  }
};rm(i, i0, r.gAll, r.jAll, r.nAll, r.gGen, r.jGen, r.nGen)

##### Regenerate x-axis in plots #####
x.axis = rep((1:(nrow(r0.gAll)/nrow(f.sce)))-1, each = nrow(f.sce))

##### Generate summary statistics #####
s0.gAll = as.data.frame(t(apply(r0.gAll, 1, summary)))
colnames(s0.gAll) = paste0("gene.", colnames(s0.gAll))
s0.jAll = as.data.frame(t(apply(r0.jAll, 1, summary)))
colnames(s0.jAll) = paste0("jump.", colnames(s0.jAll))
s0.nAll = as.data.frame(t(apply(r0.nAll, 1, summary)))
colnames(s0.nAll) = paste0("popn.", colnames(s0.nAll))
r0.plt = cbind(f.sce, x.axis, s0.gAll, s0.jAll, s0.nAll)

save(r0.plt, r0.gGen, r0.jGen, r0.nGen, sEl, sEl0, file = paste0(p.path[1],"/sum_ind.rda"))
