#!/bin/env Rscript
# author: ph-u
# script: func.r
# desc: self-defined functions for transposon ecology simulation
# in: source("func.r")
# out: NA
# arg: NA
# date: 20260623

##### Constants #####
tPn.0 = read.csv("../raw/tpn-template.csv", header = T)
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
  stop("Only four function options allowed: normal, uniform, poisson, negbin")
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
  gEne$essential = gEne$locus_tag %in% strsplit(pMs$Value[pMs$Type=="essential genes"], ";")[[1]] ## essential genes
  gEne$recombination = gEne$locus_tag %in% strsplit(pMs$Value[pMs$Type=="genes for recombination mechanism"], ";")[[1]] ## recombination mechanism
  gEne$advantage = gEne$locus_tag %in% strsplit(pMs$Value[pMs$Type=="genes for fitness advantage"], ";")[[1]] ## genes for fitness advantage
  gEnome = as.numeric(unname(gFf[gFf$type=="region", c("start", "end")]))
  gEnome = gEnome[order(gEnome)]

  ## transposon titre initiation
  tAg = "transposon population size per genome sd"
  tPn.pop = round(rNumVec(
    f = pMs$Value[pMs$Type==sub(" sd", " distribution", tAg)],
    L = pMs$Value[pMs$Type=="host organism constant population size"],
    p1 = pMs$Value[pMs$Type==sub(" sd", " mean", tAg)],
    p2 = pMs$Value[pMs$Type==tAg]),0)
  return(list(params = pMs, gene = gEne, transposon.titre = tPn.pop, genome = sum(c(gEne$length,gEne$interLength))))
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
  }else if(length(gene.var)==1){gene.var = rep(gene.var, nrow(gene.df))}
  hOst = matrix(NA,nrow = as.numeric(host.var), ncol = nrow(gene.df))
    for(i in 1:ncol(hOst)){ hOst[,i] = sample(gVar[1:gene.var[i]], size = host.var, replace = T)};rm(i)
  hOst.compact = apply(hOst, 1, function(x){paste0(x, collapse = "")})
  return(hOst.compact)
}

##### Initiate transposon pool #####
ini.transposon = function(tPn.size, scenario, template = tPn.0){
  if(length(grep(";",tPn.size))>0){tPn.size = strsplit(tPn.size, ";")[[1]]}
  tPn.size = as.numeric(tPn.size)
  tPn.wide = as.data.frame(matrix(0, nrow = length(tPn.size), ncol = ncol(template)))
  tPn.tmp = matrix(sample(LETTERS, length(tPn.size)*7, replace = T), nrow = length(tPn.size), ncol = 7)
  tPn.wide[,1] = ""
  tPn.wide[,4] = T
  tPn.wide[,5] = apply(tPn.tmp, 1, function(x){paste0(x, collapse = "")})
  tPn.wide[,6] = tPn.size
  tPn.wide[,7] = scenario$jumpRate
  tPn.wide[,8] = scenario$copyRate
  return(data.frame(ini = apply(tPn.wide, 1, function(x){paste0(x, collapse = "!")}), uniqID = tPn.wide[,5], size = tPn.size))
}

##### Rescale random number #####
reZero = function(x, new0 = 0, new1 = 1){return((x-new0)/(new1-new0))}

##### Reconstruct gene table (one host) #####
reGeneDF = function(tPn, tPn.size, gene.df){
  if(tPn != ""){
    tPn = strsplit(tPn, ";")[[1]]
    for(i in 1:length(tPn)){
      i0 = tPn.io(tPn[i])
      i1 = which(gene.df$locus_tag == substr(i0[1],2,nchar(i0[1])))
      if(substr(i0[1],1,1) == "g"){
        gene.df$length[i1] = gene.df$length[i1] + tPn.size
        if(i1<nrow(gene.df)){
          gene.df$start[(i1+1):nrow(gene.df)] = gene.df$start[(i1+1):nrow(gene.df)] + tPn.size
        }
      }else{
        gene.df$interLength[i1] = gene.df$interLength[i1] + tPn.size
        gene.df$start[i1:nrow(gene.df)] = gene.df$start[i1:nrow(gene.df)] + tPn.size
      }
      gene.df$end[i1:nrow(gene.df)] = gene.df$end[i1:nrow(gene.df)] + tPn.size
  };rm(i)}
  return(gene.df)
}

##### New host population #####
host.reproduce = function(res.pool, gene.df, transposon.size, fitness.advantage){
  # res.pool: 2 columns - $host, host genomes; $transposon, transposon notations
  ## Calculate ecological fitness deficit
  transposon.size = as.numeric(transposon.size)
  offspring.prob = rep(1, nrow(res.pool))/nrow(res.pool)
  eSsential = which(gene.df$essential)
  if(length(eSsential) > 0){ for(i in 1:length(eSsential)){ offspring.prob[grep(gene.df$locus_tag[eSsential[i]], res.pool$transposon)] = 0 };rm(i) }# dead if transposon inserted in essential genes
  aDvantage = which(gene.df$advantage) # genes with fitness advantage if transposons are inserted
  if(length(aDvantage) > 0){ for(i in 1:length(aDvantage)){ offspring.prob[grep(gene.df$locus_tag[aDvantage[i]], res.pool$transposon)] = offspring.prob[grep(gene.df$locus_tag[aDvantage[i]], res.pool$transposon)] * (1 + fitness.advantage/100) };rm(i) }
  tPn.count = lengths(strsplit(res.pool$transposon, ";"))
  tPn.count[offspring.prob==0] = 0
  offspring.prob = offspring.prob * (1 - transposon.size / sum(gene.df[,c("length", "interLength")]) * tPn.count)
  offspring.prob = offspring.prob / sum(offspring.prob)

  ## Sprouting offspring
  offspring = sample(1:nrow(res.pool), nrow(res.pool), replace = T, prob = offspring.prob)

  return(data.frame(host = res.pool$host[offspring], transposon = res.pool$transposon[offspring], familyTree = offspring))
}

##### Transposon data format conversion #####
tPn.io = function(x, ref = tPn.0){
  if(class(x)=="data.frame"){
    x$paste = apply(x,1,function(x0){paste(x0, collapse = "!")})
    return(paste(x$paste, collapse = ";"))
  }
  if(length(grep(";",x))>0){
    x = read.table(text = strsplit(x, ";")[[1]], sep = "!")
    colnames(x) = colnames(ref)
    return(x)
  }
  if(length(grep("!",x))>0){
    x = strsplit(x, "!")[[1]]
    names(x) = colnames(ref)
    return(x)
  }
  if(length(x)==ncol(ref)){return(paste(x, collapse = "!"))}
  stop("Provided ",length(x)," value(s) but ",ncol(ref)," values are needed: ",paste(colnames(ref), collapse = ", "),".")
  }

##### Modify notations according to transposon overlaps #####
tPn.x = function(tPn1, tPn2){
  t1 = tPn.io(tPn1); t2 = tPn.io(tPn2)
  if(t1[1] == t2[1]){
    t1.span = as.numeric(t1[2])*(1:as.numeric(t1[6])); t2.span = as.numeric(t2[2])*(1:as.numeric(t2[6]))
    if(as.logical(t1[4]) & any(t2.span %in% t1.span)){ # tPn1 assume always insert earlier than tPn2
      t1[4] = F
      return(tPn.io(t1))
  }}
  return(tPn1)
}

##### Revive all transposons #####
tPn.r = function(tPn){
  if(length(grep(";",tPn))>0){tPn = strsplit(tPn, ";")[[1]]}
  tPn = read.table(text = tPn, sep = "!")
  tPn[,4] = T
  return(apply(tPn, 1, function(x){paste0(x, collapse = "!")}))
}

##### Single transposon relocation #####
tPn.reloc = function(tPn, gen, gene.df){
  tPn = tPn.io(tPn)
  tPn.prob = runif(1) # set location in genome
  gene.df$cumsum = cumsum(gene.df$interLength + gene.df$length)/sum(gene.df$interLength + gene.df$length) - tPn.prob
  x = which(gene.df$cumsum >= 0)[1] # the intergenic-genic block that hosts the transposon
  if(x>1){
    tPn.loc = reZero(0, new0 = gene.df$cumsum[x-1], new1 = gene.df$cumsum[x])
  }else{
    tPn.loc = reZero(tPn.prob, new1 = gene.df$cumsum[x])
  }
  tPn.loc = ceiling(tPn.loc * (gene.df$interLength[x] + gene.df$length[x])) - gene.df$interLength[x]
  if(tPn.loc < 0){ tPn.cds = F; tPn.loc = tPn.loc + gene.df$interLength[x] }else{ tPn.cds = T }
  tPn[1:3] = c(paste0(ifelse(tPn.cds, "g", "i"), gene.df$locus_tag[x]), tPn.loc, gen)

  ## Update gene df
  gene.df$cumsum = NULL
  if(tPn.cds){
    gene.df$start[min(x+1, nrow(gene.df)):nrow(gene.df)] = gene.df$start[min(x+1, nrow(gene.df)):nrow(gene.df)] + as.numeric(tPn[6])
    gene.df$length[x] = gene.df$length[x] + as.numeric(tPn[6])
  }else{
    gene.df$start[x:nrow(gene.df)] = gene.df$start[x:nrow(gene.df)] + as.numeric(tPn[6])
    gene.df$interLength[x] = gene.df$interLength[x] + as.numeric(tPn[6])
  }
  gene.df$end[x:nrow(gene.df)] = gene.df$end[x:nrow(gene.df)] + as.numeric(tPn[6])
  colnames(gene.df)[1] = paste0(tPn.io(tPn),";",colnames(gene.df)[1])
  return(gene.df)
}

##### Hypothesized scenario modifications #####
h1.mod = function(vAl, numTpn, H1, gToxic = 1){
  x = rep(vAl,2) # prob, inherit prob -- default "fixed rate"
  if(length(grep("evolv", H1))>0){
    x = rep(abs(rnorm(1, mean = x[1], sd = .1)),2) # assume sd = 0.1
  }else if(length(grep("charlesworth", H1))>0){
    x[1] = x[1]/numTpn
  }
  x[1] = x[1]*gToxic # genotoxic effect not inheritable
  return(x)
}

sCene.mod = function(pRobs, tPn.bg, sCene, pAram, gToxic = F){
  tPn.bg = length(strsplit(tPn.bg, ";")[[1]])
  gTx = runif(2) < pAram[1:2]
  x = c(
    h1.mod(pRobs[1], tPn.bg, sCene[2], ifelse(gToxic, ifelse(gTx[1], 0, ifelse(gTx[2], pAram[3], 1)),1)), # jump
    h1.mod(pRobs[2], tPn.bg, sCene[4], 1) # copy is not affected by genotoxic effect
  )
  return(x[c(1,3,2,4)])
}

##### Get transposons from gene.df #####
tPn.get = function(gene.df, transposon = T){
  if(transposon){
    return(paste(rev(rev(strsplit(colnames(gene.df)[1], ";")[[1]])[-1]), collapse = ";"))
  }else{
    return(rev(strsplit(colnames(gene.df)[1], ";")[[1]])[1])
  }
}

##### Single transposon action: jump and/or copy or neither? #####
tPn.act = function(tPn, gen, gene.df, scenario, pAram, gToxic = F){ # gene.df must have all transposons attached
  x = tPn.io(tPn); x0 = tPn
  x.probs = sCene.mod(as.numeric(x[c("jump", "copy")]), tPn.get(gene.df), scenario, pAram, gToxic)
  x[c("jump", "copy")] = x.probs[3:4]
  a = c(rNumVec(f = "uniform", L = 2, p1 = 0, p2 = 1) < x.probs[1:2], max(1, rpois(1, rnorm(1, 1, .01)))) # jump?, copy?, num copies
  a[1] = ifelse(gen > 0,a[1],1);a[2] = ifelse(gen > 0,a[2],0)
  colnames(gene.df)[1] = tPn.get(gene.df, F)
  if(as.logical(x["valid"])){
    if(a[1]>0){ # jump
      gene.df = tPn.reloc(tPn, gen, gene.df)
      x0 = tPn.get(gene.df)
      colnames(gene.df)[1] = tPn.get(gene.df, F)
    }
    if(a[2]>0){ # copy
      x0 = unique(c(tPn,x0))
      for(i in 1:a[3]){
      gene.df = tPn.reloc(tPn, gen, gene.df)
      tTpn = c(tPn, tPn.get(gene.df))
      tCheck = abs(unlist(apply(tPn.io(paste(tTpn, collapse = ";")), 1, function(x){which(gene.df[,1]==substr(x[1], 2, nchar(x[1])))})) - nrow(gene.df)/2)
      if(scenario$copyDir == "terminus"){ # assume transposon is not quick enough to copy into the adjacent intragenic-gene pair
        cpValid = tCheck[2] < tCheck[1]
      }else if(scenario$copyDir == "origin"){
        cpValid = tCheck[2] > tCheck[1]
      }else{ cpValid = T }
      if(cpValid){
        x0 = c(x0, tPn.get(gene.df))
      }
      colnames(gene.df)[1] = tPn.get(gene.df, F)
      }}
    colnames(gene.df)[1] = paste0(c(x0, colnames(gene.df)[1]), collapse = ";")
  }
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
g.Recom = function(res.pool, gene.df, recomRate, hypothesis = "switch"){
  # res.pool: 2 columns - $host, host genomes; $transposon, transposon notations
  if(hypothesis == "homeostatic"){
    numGenes = floor(lengths(strsplit(res.pool$transposon, ";")) * as.numeric(recomRate) * nrow(gene.df)) # gene recombination rate has a linear increase according to the number of transposons in the genome (homeostatic)
  }else{
    numGenes = rep(floor(as.numeric(recomRate) * nrow(gene.df)), nrow(res.pool)) # switch (default)
  }

  ## If any gene in recombination mechanism is hit by transposon, no recombination
  if(sum(gene.df$recombination)>0){ for(i in which(gene.df$recombination)){
    numGenes[grep(gene.df$locus_tag[i], res.pool$transposon)] = 0
  };rm(i)}
  numGenes[numGenes > (nrow(gene.df)-1)] = nrow(gene.df)-1 # random selection of genes can't exceed total number of genes in genome

  for(i in 1:nrow(res.pool)){ if(numGenes[i] > 0){
    x.recom = data.frame(
      gene = sample(c(1:nrow(gene.df))[-i], numGenes[i], replace = F), # which gene being recombined
      host = sample(c(1:nrow(res.pool))[-i], numGenes[i], replace = T)) # which host is the source
    for(i0 in 1:nrow(x.recom)){
      x.receipt = gene.recom(h1.G = res.pool$host[i], h1.t = res.pool$transposon[i], h2.G = res.pool$host[x.recom$host[i0]], h2.t = res.pool$transposon[x.recom$host[i0]], g2to1 = x.recom$gene[i0], locusTags = gene.df$locus_tag)
      res.pool$host[i] = x.receipt[1]
      res.pool$transposon[i] = sub(";$","",sub("^;","",x.receipt[2]))
    };rm(i0, x.receipt, x.recom)
  }};rm(i)
  return(res.pool)
}
