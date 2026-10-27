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
src.url = "https://raw.githubusercontent.com/ph-u/snippet_ph-u/refs/heads/master/"
source(paste0(src.url, "colour.r"))

##### f: Decode transposon annotation (from func.r) #####
#tPn.io = function(x, ref = f.tpn){
#  if(class(x)=="data.frame"){
#    x$paste = apply(x,1,function(x0){paste(x0, collapse = "!")})
#    return(gsub(" ", "", paste(x$paste, collapse = ";")))
#  }
#  if(length(grep(";",x))>0){
#    x = read.table(text = strsplit(x, ";")[[1]], sep = "!")
#    colnames(x) = colnames(ref)
#    return(x)
#  }
#  if(length(grep("!",x))>0){
#    x = strsplit(x, "!")[[1]]
#    names(x) = colnames(ref)
#    return(x)
#  }
#  if(length(x)==ncol(ref)){return(gsub(" ", "", paste(x, collapse = "!")))}
#  stop("Provided ",length(x)," value(s) but ",ncol(ref)," values are needed: ",paste(colnames(ref), collapse = ", "),".")
#}

##### f: Calculate percentage contribution of source #####
# famRoot = function(p1 = "1!0.5;2!0.3;3!0.2", p2 = "1!0.5;2!0.3;3!0.2"){
#   p1 = read.table(text = strsplit(p1, ";")[[1]], sep = "!")
#   p2 = read.table(text = strsplit(p2, ";")[[1]], sep = "!")
#   p1[,2] = (p1[,2] + p2[,2])/2
#   return(paste(p1[,1], p1[,2], sep = "!", collapse = ";"))
# }

##### f: Reconstruct family tree #####
fAmily = function(tree.df, state = "haploid"){
  nGen = nrow(tree.df); nPop = ncol(tree.df)
  par = array(0, c(2, nPop, nGen))
  for(i in seq_len(nGen)){par[,,i] = matrix(as.integer(unlist(strsplit(unlist(tree.df[i,], use.names=FALSE), ";"), use.names=FALSE)), nrow = 2)}
  src = sort(unique(as.vector(par[,,1]))); K = length(src)
  out = matrix(0, nGen, K)
  if(state=="haploid"){
    P = matrix(0, K, nPop); j = seq_len(nPop)
    for(k in 1:2){ idx = cbind(match(par[k,,1], src), j); P[idx] = P[idx] + .5 }
    out[1,] = rowSums(P)
    if(nGen > 1){ for(i in 2:nGen){
      P = (P[, par[1,,i] - ifelse(par[1,,i] > nPop, nPop, 0), drop=F] + P[, par[2,,i] - ifelse(par[2,,i] > nPop, nPop, 0), drop=F])/2
      out[i,] = rowSums(P)
    }}
    as.data.frame(out)
  }else{
    out[1,] = tabulate(P <- match(as.numeric(t(par[,,1])), src), nbins = K)
    if(nGen > 1){for(i in 2:nGen){ out[i,] = tabulate(P <- P[as.numeric(t(par[,,i]))], nbins = K) }}
    as.data.frame(out/dim(par)[1])
  }
}
