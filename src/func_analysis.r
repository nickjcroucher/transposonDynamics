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
famRoot = function(p1 = "1!0.5;2!0.3;3!0.2", p2 = "1!0.5;2!0.3;3!0.2"){
  p1 = read.table(text = strsplit(p1, ";")[[1]], sep = "!")
  p2 = read.table(text = strsplit(p2, ";")[[1]], sep = "!")
  p1[,2] = (p1[,2] + p2[,2])/2
  return(paste(p1[,1], p1[,2], sep = "!", collapse = ";"))
}

##### f: Reconstruct family tree #####
fAmily = function(tree.df){
  nGen = nrow(tree.df); nPop = ncol(tree.df)
  par = array(0, c(2, nPop, nGen))
  for(i in seq_len(nGen)){par[,,i] = matrix(as.integer(unlist(strsplit(unlist(tree.df[i,], use.names=FALSE), ";"), use.names=FALSE)), nrow = 2)}
  src = sort(unique(as.vector(par[,,1]))); K = length(src)
  P = matrix(0, K, nPop); j = seq_len(nPop)
  for(k in 1:2){ idx = cbind(match(par[k,,1], src), j); P[idx] = P[idx] + .5 }
  out = matrix(0, nGen, K); out[1,] = rowSums(P)
  for(i in 2:nGen){
    P = (P[, par[1,,i], drop=F] + P[, par[2,,i], drop=F])/2
    out[i,] = rowSums(P)
  }
  as.data.frame(out)
}

##### f: Get transposon generations #####
tpn.gen = function(x, ref = tPn.0){
  nG = nrow(x); nP = ncol(x); nF = ncol(ref)
  iGen = which(colnames(ref) == "generation")
  cnt  = matrix(0, nG, nP)
  tgen = matrix(0, nG, nG)
  dist = matrix(0, nG, 2, dimnames = list(NULL, c("g","i")))
  for(i in seq_len(nG)){
    cells = as.character(unlist(x[i,], use.names = F))
    keep  = nzchar(cells)
    if(!any(keep)) next
    tok.l   = strsplit(cells, ";", fixed = T)
    cnt[i,] = lengths(tok.l) * keep
    tok = unlist(tok.l, use.names = F); tok = tok[nzchar(tok)]
    if(!length(tok)) next
    sp = strsplit(tok, "!", fixed = T)
    if(any(lengths(sp) != nF)){ stop("malformed transposon record at generation ", i-1) }
    fld = matrix(unlist(sp, use.names = F), nrow = nF)
    tgen[i,] = tabulate(as.integer(fld[iGen,]) + 1, nbins = nG)
    pfx = substr(tok, 1, 1)
    dist[i,] = c(sum(pfx == "g"), sum(pfx == "i"))
  }
  return(list(count = as.data.frame(cnt), generation = as.data.frame(tgen), distribution = as.data.frame(dist)))
}

##### f: Default loading simulation result file #####
simLoad = function(x){
  load(x)
  tpn = tpn.gen(rec.transposon)
  return(list(host = rec.host, transposon = rec.transposon, tpn.count = tpn$count, tpn.generation = tpn$generation, tpn.distribution = tpn$distribution, genealogy = fAmily(rec.offspring)/ncol(rec.offspring)))
}
