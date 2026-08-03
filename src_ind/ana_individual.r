#!/bin/env Rscript
# author: ph-u
# script: ana_individual.r
# desc: analysis from point-of-view of one genome
# in: Rscript ana_individual.r
# out: NA
# arg: 0
# date: 20260716

argv = as.numeric(commandArgs(T))
if(length(argv) < 1){argv = 1}

##### set env #####
if(Sys.info()["sysname"]=="Darwin"){
  source("../src/p_path.r")
}else{
  p.file = list.files("..", full.names = T, recursive = T)
}
source("../src/func_analysis.r")
a.rDa = p.file[grep("tPn-", p.file)]
a.rDa = a.rDa[grep("rda$", a.rDa)]
if(Sys.info()["sysname"]=="Darwin"){ a.rDa = a.rDa[grep("_pj02", a.rDa)] }
f.rda = read.table(text = sub("[.]", "_", read.table(text = a.rDa, sep = "-")[,3]), sep = "_")[,1:2]
colnames(f.rda) = c("rep", "sce")
f.rda$file = a.rDa;rm(a.rDa)
f.sce = read.csv(p.file[grep("scenario", p.file)[1]], header = T)
f.hos = read.csv(p.file[grep("template-host", p.file)[1]], header = T)
f.tpn = read.csv(p.file[grep("template-tpn", p.file)[1]], header = T)
#f.sce$paste = apply(f.sce, 1, function(x){paste0(x, collapse = "")})
#f.rda = cbind(f.rda, f.sce[f.rda$sce,])
#f.rda$hypothesis = apply(f.rda, 1, function(x){paste0(x[6:7], collapse = "")})

##### Set data collection structure templates #####
#options(warn = 2)
#set.seed(1234)
s1 = simLoad(f.rda$file[1])

r.res0 = as.data.frame(matrix(nrow = nrow(rec.transposon)*nrow(f.sce), ncol = max(f.rda$rep))) # overall result dataframe
r.res = cbind(f.sce, r.res0[,1:2]); r.res[,(-1:0)+ncol(r.res)] = NULL
r.res$generation = rep((1:nrow(rec.transposon))-1, each = nrow(f.sce))
#r.res$hypothesis = paste0(r.res$jump,r.res$recom)

cat(date(), ": collecting individualized resolution data\n")

## Pick sample replicates
sEl = c();sEl0 = unique(r.res$hypothesis)
r.sel = vector(mode = "list", length = length(sEl0))
for(i in 1:length(sEl0)){
  sEl = c(sEl, sample(which(f.rda$hypothesis == sEl0[i]), 1))
  r.sel[[i]] = as.data.frame(matrix(nrow = nrow(rec.transposon), ncol = ncol(rec.transposon)))
};rm(i)

##### Data-collection dataframes #####
r.nAll = r.res0; r.nGen = r.sel # titre
r.jAll = r.res0; r.jGen = r.sel # jumping generations
r.gAll = r.res0; r.gGen = r.sel # gene/intragenic distributions

##### Data collection #####
for(i in (argv[1]-1)*5+(1:5)){ #1:nrow(f.rda)){
  cat(date(), ":", i, "/", nrow(f.rda), "(", round(i/nrow(f.rda)*100, 2), "% )      \n")
  ok = tryCatch({
    load(f.rda$file[i])
    T
  }, error = function(s){cat(date(),":",f.rda$file[i],"data corruption error (",i,")\n")}, warning = function(s){cat(date(),":",f.rda$file[i],"data corruption warning (",i,")\n")})
  if(length(ok)>0 && ok){
    if(i>1){ load(f.rda$file[i]) }
    for(i0 in 1:nrow(rec.transposon)){
#      cat(date(), ":", i, "/", nrow(f.rda), "(", round(i/nrow(f.rda)*100, 2), "% );", i0, "/", nrow(rec.transposon), "(", round(i0/nrow(rec.transposon)*100, 2), "% )      \n")
      r.tmp = strsplit(as.character(rec.transposon[i0,]), ";")
      r.loc = c(which(r.res$paste == f.rda$paste[i] & r.res$generation == (i0-1)),paste0("V",f.rda$rep[i]))

## Overall individualized distributions
      r.nAll[as.numeric(r.loc[1]), r.loc[2]] = mean(lengths(r.tmp)) # mean titre

      r.gG0 = rep(NA, length(r.tmp))
      for(i1 in 1:length(r.tmp)){
        if(length(grep("r.jG0", objects()))>0 && i1==1){rm(r.jG0)}
        if(length(r.tmp[[i1]])>0){
          r.00 = read.table(text = r.tmp[[i1]], sep = "!")
## jumping generations
          r.jG1 = as.data.frame(table(r.00[,3]))
          colnames(r.jG1)[2] = paste0(colnames(r.jG1)[2],i1)
          if(i1>1){r.jG0 = merge(r.jG0, r.jG1, all = T)}else{r.jG0 = r.jG1}

## gene distributions
          r.gG0[i1] = sum(substr(r.00[,1],1,1)=="g")/nrow(r.00)
        }else if(i1==1){
          r.jG0 = as.data.frame(table(0))
          colnames(r.jG0)[2] = paste0(colnames(r.jG0)[2],i1)
	  r.jG0[1,2] = 0
        }else{
          r.jG0[,(ncol(r.jG0)+1)] = NA
          colnames(r.jG0)[ncol(r.jG0)] = paste0("Freq",i1)
      }};rm(i1)
      r.jG0[is.na(r.jG0)] = 0
      r.jG0$sum = rowSums(r.jG0[,-1])
      r.jG0$sum = r.jG0$sum/sum(r.jG0$sum)
      r.jAll[as.numeric(r.loc[1]), r.loc[2]] = sum(as.numeric(r.jG0[,1])*r.jG0$sum)# mean jumping generations
      r.gAll[as.numeric(r.loc[1]), r.loc[2]] = mean(r.gG0) # mean gene distributions

## Collection for selected replicates
      if(i %in% sEl){
        cat(date(),": capturing ",i,"(tag num) simulation\n")
        i.loc = which(sEl == i)
        r.nGen[[i.loc]][i0,] = lengths(r.tmp)
	r.jG0$sum = NULL
	for(i1 in 2:ncol(r.jG0)){r.jG0[,i1] = as.numeric(r.jG0[,1])*r.jG0[,i1]/sum(r.jG0[,i1])}
	r.jGen[[i.loc]][i0,] = colSums(r.jG0[,-1])
        r.gGen[[i.loc]][i0,] = r.gG0
      }
    }
}else{cat(date(),":",f.rda$file[i],"data corruption\n")}
};cat(date(), ": done\n")
rm(i,i0, rec.host, rec.offspring, rec.transposon)

save(r.nAll, r.nGen, r.jAll, r.jGen, r.gAll, r.gGen, sEl, sEl0, f.rda, f.sce, file = paste0("../data/ana_ind-",argv[1],".rda"))
