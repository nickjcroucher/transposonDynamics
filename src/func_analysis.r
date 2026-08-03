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
  cat(date(), ": mapping family tree\n")
  n.var = as.numeric(unique(strsplit(paste0(tree.df[1,], collapse = ";"), ";")[[1]]))
  n.var.df = data.frame(src = n.var[order(n.var)], prop = 0)
  for(i in 1:length(n.var)){
    n.0 = n.var.df; n.0[i,2] = 1
    n.var[i] = paste(n.0[,1], n.0[,2], sep = "!", collapse = ";")
  };rm(i)
  t.df = as.data.frame(matrix(0, nrow = nrow(tree.df), ncol = length(n.var)))
  for(i in 1:nrow(tree.df)){
#cat(date(), ":", i, "/", nrow(tree.df), "(", round(i/nrow(tree.df)*100,2), "% )     \r")
    for(i0 in 1:ncol(tree.df)){
      i1 = as.numeric(strsplit(tree.df[i,i0], ";")[[1]])
      if(i > 1){
        tree.df[i,i0] = famRoot(p1 = tree.df[i-1,i1[1]], p2 = tree.df[i-1,i1[2]])
      }else{
        tree.df[i,i0] = famRoot(p1 = n.var[i1[1]], p2 = n.var[i1[2]])
      }
      if(i0 > 1){
        i2[,2] = i2[,2] + read.table(text = strsplit(tree.df[i,i0], ";")[[1]], sep = "!")[,2]
      }else{
        i2 = read.table(text = strsplit(tree.df[i,i0], ";")[[1]], sep = "!")
      }
    };rm(i0, i1)
    t.df[i,] = i2[,2]
  };rm(i2);cat("\n")
  t.df[is.na(t.df)] = 0
  cat(date(), ": mapping done\n")
  return(t.df)
}

##### f: Get transposon generations #####
tpn.gen = function(x){
  t.df = as.data.frame(matrix(0, nrow = nrow(x), ncol = nrow(x)))
  g.df = as.data.frame(matrix(0, nrow = nrow(x), ncol = 2))
  colnames(g.df) = c("g", "i")
  n.df = as.data.frame(matrix(0, nrow = nrow(x), ncol = ncol(x)))
  for(i in 1:nrow(x)){ if(sum(as.character(x[i,])!="")>0){
    n.df[i,] = lengths(strsplit(as.character(x[i,]), ";"))
    i0 = tPn.io(paste0(x[i,], collapse=";"))
    i1 = table(i0$generation)
    t.df[i,as.numeric(names(i1))+1] = i1
    i1 = table(substr(i0$gene,1,1))
    g.df[i,names(i1)] = i1
  }};rm(i, i0, i1)
  return(list(count = n.df, generation = t.df, distribution = g.df))
}

##### f: Default loading simulation result file #####
simLoad = function(x){
  load(x)
  tpn = tpn.gen(rec.transposon)
  return(list(host = rec.host, transposon = rec.transposon, tpn.count = tpn$count, tpn.generation = tpn$generation, tpn.distribution = tpn$distribution, genealogy = fAmily(rec.offspring)/ncol(rec.offspring)))
}
