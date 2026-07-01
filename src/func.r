#!/bin/env Rscript
# author: ph-u
# script: func.r
# desc: self-defined functions for transposon ecology simulation
# in: source("func.r")
# out: NA
# arg: NA
# date: 20260623

##### Constants #####
# write.csv(data.frame(transposon = rep(10^-(2:7), each = 6), gene = rep(10^-(2:7), 6)), "../raw/scenario.csv", row.names = F, quote = F)
gVar = c(letters,LETTERS,0:9) # notations for genome variation

##### Reorganize GFF file information #####
gffClean = function(gFf){
  colnames(gFf) = c("seqname","source","type","start","end","score","strand","phase","attributes")
  x = strsplit(gFf$attributes, ";")
  x0 = unlist(strsplit(sub("=","=@",unlist(x)), "="))
  cNam = unique(x0[grep("@", x0, invert = T)])
  aTt = as.data.frame(matrix(nrow = length(x), ncol = length(cNam)))
  colnames(aTt) = cNam
  for(i0 in 1:length(x)){for(i1 in 1:length(x[[i0]])){
    i2 = strsplit(x[[i0]][i1], "=")[[1]]
    aTt[i0,colnames(aTt)==i2[1]] = i2[2]
  }};rm(i0,i1)
  gFf = cbind(gFf[,-ncol(gFf)],aTt)
  return(gFf)
}

##### Random number vector generator #####
rNumVec = function(f="", L=1, p1=1, p2=0){
  L = as.numeric(L); p1 = as.numeric(p1); p2 = as.numeric(p2)
  if(f=="normal"){ return(rnorm(n = L, mean = p1, sd = p2))
  }else if(f=="uniform"){ return(runif(n = L, min = max(0, p1-p2), max = p1+p2))
  }else if(f=="poisson"){ return(1/(rpois(n = L, lambda = p1) + 1))
  }else if(f=="negbin"){ return(1/(rnbinom(n = L, prob = p1, size = p2) + 1))}
  stop("Only three function options allowed: normal, uniform, poisson")
}

##### Format input parameters, gene table, genome stretch, and transposon titres in a host population #####
inParams = function(pArams = "../raw/input.csv"){
  pMs = read.csv(pArams, header = T)
  gFf = gffClean(read.table(paste0("../data/",strsplit(pMs$Value[pMs$Type=="ref genome"], "[.]")[[1]][1],".gff"), sep = "\t", header = F, quote = ""))
  gEne = unique(gFf[!is.na(gFf$product),c("locus_tag","start","end","product")])
  gEne$length = gEne$end - gEne$start + 1
  gEne$interLength = c(gEne$start[1]-1+gFf[gFf$type=="region", "end"]-gEne$end[nrow(gEne)], gEne$start[-1] - gEne$end[-nrow(gEne)])
  x = which(gEne$interLength<0)
  gEne$length[x-1] = gEne$length[x-1] + gEne$interLength[x]
  gEne$interLength[x] = 0
  gEne$essential = F ## + essential genes
  gEne$recombination = F ## + recombination mechanism
  gEnome = as.numeric(unname(gFf[gFf$type=="region", c("start", "end")]))
  gEnome = gEnome[order(gEnome)]

  ## transposon titre initiation
  tAg = "transposon population size per genome sd"
  tPn.pop = round(1/rNumVec(
    f = pMs$Value[pMs$Type==sub(" sd", " distribution", tAg)],
    L = pMs$Value[pMs$Type=="host organism constant population size"],
    p1 = pMs$Value[pMs$Type==sub(" sd", " mean", tAg)],
    p2 = pMs$Value[pMs$Type==tAg]),0)
  return(list(params = pMs, gene = gEne, transposon.titre = tPn.pop))
}

##### Initiate genome pool #####
ini.host = function(host.var, gene.df, gene.var){
  if(length(grep(";", gene.var))>0){gene.var = strsplit(gene.var, ";")[[1]]}
  gene.var = as.numeric(gene.var)
  if(!(length(gene.var) %in% c(1,nrow(gene.df)))){
    if(length(gene.var) > nrow(gene.df)){
      warning(paste0("Too many notations, extra ones are ignored. You provided ",length(gene.var), " gene variation notations, your gff file only has ",nrow(gene.df), " genes."))
      gene.var = gene.var[1:nrow(gene.df)]
    }else{
      stop(paste0("Either one variation or specify variations for each gene in the genome. You only provided ",length(gene.var), " gene variation notations, your gff file has ",nrow(gene.df), " genes."))
    }
  }
  hOst = matrix(NA,nrow = as.numeric(host.var), ncol = nrow(gene.df))
    for(i in 1:ncol(hOst)){ hOst[,i] = sample(gVar[1:gene.var[i]], size = host.var, replace = T)};rm(i)
  hOst.compact = apply(hOst, 1, function(x){paste0(x, collapse = "")})
  return(hOst.compact)
}

##### Initiate transposon pool #####
ini.transposon = function(tPn.size){
  if(length(grep(";",tPn.size))>0){tPn.size = strsplit(tPn.size, ";")[[1]]}
  tPn.size = as.numeric(tPn.size)
  tPn.wide = as.data.frame(matrix(0, nrow = length(tPn.size), ncol = 5))
  tPn.tmp = matrix(sample(LETTERS, length(tPn.size)*7, replace = T), nrow = length(tPn.size), ncol = 7)
  tPn.wide[,ncol(tPn.wide)] = apply(tPn.tmp, 1, function(x){paste0(x, collapse = "")})
  tPn.wide[,4] = "T"
  tPn.wide[,1] = ""
  return(data.frame(ini = apply(tPn.wide, 1, function(x){paste0(x, collapse = "!")}), uniqID = tPn.wide[,ncol(tPn.wide)], size = tPn.size))
}

##### Rescale random number #####
reZero = function(x, new0 = 0, new1 = 1){return((x-new0)/(new1-new0))}

##### Reconstruct gene table (one host) #####
reGeneDF = function(tPn, tPn.size, gene.df){
  if(tPn != ""){
    tPn = strsplit(tPn, ";")[[1]]
    for(i in 1:length(tPn)){
      i0 = tPn.io(chain = tPn[i], encrypt = F)
      i1 = which(gene.df$locus_tag == substr(i0[1],2,nchar(i0[1])))
      if(substr(i0[1],1,1) == "g"){
        gene.df$length[i1] = gene.df$length[i1] + tPn.size
        gene.df$start[(i1+1):nrow(gene.df)] = gene.df$start[(i1+1):nrow(gene.df)] + tPn.size
      }else{
        gene.df$interLength[i1] = gene.df$interLength[i1] + tPn.size
        gene.df$start[i1:nrow(gene.df)] = gene.df$start[i1:nrow(gene.df)] + tPn.size
      }
      gene.df$end[i1:nrow(gene.df)] = gene.df$end[i1:nrow(gene.df)] + tPn.size
  };rm(i)}
  return(gene.df)
}

##### New host population #####
host.reproduce = function(res.pool, gene.df, transposon.size){
  # res.pool: 2 columns - $host, host genomes; $transposon, transposon notations
  ## Calculate ecological fitness deficit
  transposon.size = as.numeric(transposon.size)
  offspring.prob = rep(1, nrow(res.pool))/nrow(res.pool)
  eSsential = which(gene.df$essential)
  if(length(eSsential) > 0){ # dead if transposon inserted in essential genes
    for(i in 1:length(eSsential)){
      dEad = grep(gene.df$locus_tag[eSsential[i]], res.pool$transposon)
      extra.prob = sum(offspring.prob[dEad])
      offspring.prob[dEad] = 0
      offspring.prob[-dEad] = offspring.prob[-dEad] + extra.prob/length(offspring.prob[-dEad])
    };rm(i,dEad,extra.prob)
  }
  tPn.count = lengths(strsplit(res.pool$transposon, ";"))
  tPn.count[offspring.prob==0] = 0
  offspring.prob = offspring.prob * (1 - transposon.size / sum(gene.df[,c("length", "interLength")]) * tPn.count)
  offspring.prob = offspring.prob / sum(offspring.prob)

  ## Sprouting offspring
  offspring = sample(1:nrow(res.pool), nrow(res.pool), replace = T, prob = offspring.prob)

  return(data.frame(host = res.pool$host[offspring], transposon = res.pool$transposon[offspring], familyTree = offspring))
}

##### Transposon data format conversion #####
tPn.io = function(gene = "", location = "", generation = "", valid = "", uniqID = "", chain = "", encrypt = T){
  if(encrypt==T){
    x = paste(c(gene, location, generation, valid, uniqID), collapse = "!")
  }else{
    x = strsplit(chain, "!")[[1]]
  };return(x)
}

##### Modify notations according to transposon overlaps #####
tPn.x = function(tPn1, tPn1.len, tPn2, tPn2.len){
  t1 = tPn.io(chain = tPn1, encrypt = F); t2 = tPn.io(chain = tPn2, encrypt = F)
  t1.span = as.numeric(t1[2]):(as.numeric(t1[2])+as.numeric(tPn1.len)); t2.span = as.numeric(t2[2]):(as.numeric(t2[2])+as.numeric(tPn2.len))
  if(t1[4]=="T" && any(t2.span %in% t1.span)){ # tPn1 assume always insert earlier than tPn2
    return(tPn.io(gene = t1[1], location = t1[2], generation = t1[3], valid = "F", uniqID = t1[5], encrypt = T))
  }else{
    return(tPn1)
  }
}

##### Revive all transposons #####
tPn.r = function(tPn){
  if(length(grep(";",tPn))>0){tPn = strsplit(tPn, ";")[[1]]}
  tPn = read.table(text = tPn, sep = "!", colClasses = "character")
  tPn[,4] = "T"
  tPn.ind = apply(tPn, 1, function(x){paste0(x, collapse = "!")})
  return(paste0(tPn.ind, collapse = ";"))
}

##### Single transposon jump #####
tPn.jump = function(tPn.tag, tPn.prob, tPn.size, gene.df, tPn.move = 0, jumProb = 0, generation = 1){
  if(length(grep(".", colnames(gene.df)[1], fixed = T))>0){colnames(gene.df)[1] = strsplit(colnames(gene.df)[1], "[.]")[[1]][2]} # clear previous tPn.tag
  if(tPn.move < jumProb){ # Determine location of jump
    tPn.prob = reZero(tPn.prob, new1 = jumProb)
    gene.df$cumsum = cumsum(gene.df$interLength + gene.df$length)/sum(gene.df$interLength + gene.df$length) - tPn.prob
    x = which(gene.df$cumsum > 0)[1] # the intergenic-genic block that hosts the transposon
    if(x>1){
      tPn.loc = reZero(0, new0 = gene.df$cumsum[x-1], new1 = gene.df$cumsum[x])
    }else{
      tPn.loc = reZero(tPn.prob, new1 = gene.df$cumsum[x])
    }
    tPn.loc = ceiling(tPn.loc * (gene.df$interLength[x] + gene.df$length[x])) - gene.df$interLength[x]
    if(tPn.loc < 0){ tPn.cds = F; tPn.loc = tPn.loc + gene.df$interLength[x] }else{ tPn.cds = T }

    tPn.old = tPn.io(chain = tPn.tag, encrypt = F)
    tPn.tag = tPn.io(gene = paste0(ifelse(tPn.cds, "g", "i"), gene.df$locus_tag[x]), location = tPn.loc, generation = generation, valid = "T", uniqID = tPn.old[5])

    ## Update gene df
    gene.df$cumsum = NULL
    if(tPn.cds){
      gene.df$start[min(x+1, nrow(gene.df)):nrow(gene.df)] = gene.df$start[min(x+1, nrow(gene.df)):nrow(gene.df)] + tPn.size
    }else{
      gene.df$start[x:nrow(gene.df)] = gene.df$start[x:nrow(gene.df)] + tPn.size
    }
    gene.df$end[x:nrow(gene.df)] = gene.df$end[x:nrow(gene.df)] + tPn.size
  }
  colnames(gene.df)[1] = paste0(tPn.tag,".",colnames(gene.df)[1])
  return(gene.df)
}

##### Single gene recombination #####
gene.recom = function(h1.G, h1.t, h2.G, h2.t, g2to1, locusTags){
  h1.G = strsplit(h1.G,"")[[1]]
  h1.G[g2to1] = substr(h2.G,g2to1,g2to1)
  h1.G = paste0(h1.G, collapse = "")
  if(length(grep(locusTags[g2to1], h1.t)) > 0){
    h1.t = strsplit(h1.t, ";")[[1]]
    h1.t = h1.t[-grep(locusTags[g2to1], h1.t)]
    h1.t = paste0(h1.t, collapse = ";")
  }
  if(length(grep(locusTags[g2to1], h2.t)) > 0){
    h2.t = strsplit(h2.t, ";")[[1]]
    h1.t = paste0(c(h1.t, h2.t[grep(locusTags[g2to1], h2.t)]), collapse = ";")
  }
  return(c(genome = h1.G, transposon = h1.t))
}

##### Gene recombination in host population, assume ascending order as recipients #####
g.Recom = function(res.pool, gene.df, recomRate){
  # res.pool: 2 columns - $host, host genomes; $transposon, transposon notations
  numGenes = floor(lengths(strsplit(res.pool$transposon, ";")) * as.numeric(recomRate) * nrow(gene.df)) # gene recombination rate has a linear increase according to the number of transposons in the genome
  numGenes[numGenes > (nrow(gene.df)-1)] = nrow(gene.df)-1
  for(i in 1:nrow(res.pool)){ if(numGenes[i] > 0){
    x.recom = data.frame(
      gene = sample(c(1:nrow(gene.df))[-i], numGenes[i], replace = F), # which gene being recombined
      host = sample(c(1:nrow(res.pool))[-i], numGenes[i], replace = T)) # which host is the source
    for(i0 in 1:nrow(x.recom)){
      x.receipt = gene.recom(h1.G = res.pool$host[i], h1.t = res.pool$transposon[i], h2.G = res.pool$host[x.recom$host[i0]], h2.t = res.pool$transposon[x.recom$host[i0]], g2to1 = x.recom$gene[i0], locusTags = gene.df$locus_tag)
      res.pool$host[i] = x.receipt[1]
      res.pool$transposon[i] = sub(";$","",sub("^;","",x.receipt[2]))
    };rm(i0, x.receipt)
  }};rm(i, x.recom)
  return(res.pool)
}
