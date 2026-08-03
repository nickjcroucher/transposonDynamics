#!/bin/env Rscript
# author: ph-u
# script: func_analysis.r
# desc: self-defined functions for analyses on the transposon ecology simulations
# in: source("func_analysis.r")
# out: NA
# arg: NA
# date: 20260703

#for(i in c("vegan", "fossil", "gganimate", "ggplot2", "ggtern")){suppressPackageStartupMessages(library(i, character.only = T))};rm(i)
#source("func.r")

##### Colour #####
cBp = cBl = cBl0 = c();for(i in c("Okabe-Ito", "Alphabet", "Polychrome 36", "Classic Tableau")){
  cBp = c(cBp,palette.colors(palette = i, alpha=1, recycle = F))
  cBl = c(cBl,palette.colors(palette = i, alpha=.7, recycle = F))
  cBl0 = c(cBl0,palette.colors(palette = i, alpha=.3, recycle = F))
};rm(i)

##### f: Reconstruct family tree #####
fAmily = function(tree.df){
  cat(date(), ": mapping family tree\n")
  n.var = as.numeric(unique(strsplit(paste0(tree.df[1,], collapse = ";"), ";")[[1]]))
  n.var = n.var[order(n.var)]
  t.df = as.data.frame(matrix(0, nrow = nrow(tree.df), ncol = length(n.var)))
  for(i in 1:nrow(tree.df)){
    cat(date(), ":", i, "/", nrow(tree.df), "(", round(i/nrow(tree.df)*100,2), "% )     \r")
    i0 = as.data.frame(table(strsplit(paste0(tree.df[i,], collapse = ";"), ";")[[1]]))
    i0[,1] = as.character(i0[,1])
    if(i0>1){
      i0[,3] = as.character(tree.df[i-1,as.numeric(i0[,1])])
    }
    t.df[i,] = i0$Freq[match(n.var, i0$Var1)]/sum(i0$Freq)
  };rm(i);cat("\n")
  t.df[is.na(t.df)] = 0
  cat(date(), ": mapping done\n")
  return(t.df)}

##### f: Default loading simulation result file #####
simLoad = function(x){
  load(x)
  fAm = fAmily(rec.offspring)
  return(list())
}
