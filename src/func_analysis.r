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

##### Reconstruct family tree #####
fAmily = function(tree.df){
  for(i in 2:nrow(tree.df)){
    tree.df[i,] = as.numeric(tree.df[i-1,as.numeric(tree.df[i,])])
  };rm(i)
  return(tree.df)
}
