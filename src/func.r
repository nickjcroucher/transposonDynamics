#!/bin/env Rscript
# author: ph-u
# script: func.r
# desc: self-defined functions for transposon ecology simulation
# in: source("func.r")
# out: NA
# arg: NA
# date: 20260623

##### Constants #####
tPn.0 = read.csv("../raw/template-tpn.csv", header = T)
host.0 = read.csv("../raw/template-host.csv", header = T)

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
  }else if(f=="uniform"){ return(runif(n = L, min = max(0, min(1, p1-p2)), max = min(1, p1+p2)))
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
  gEne = rbind(gEne,gEne)
  gEne$locus_tag = paste0(rep(0:1, each = nrow(gEne)/2), gEne$locus_tag)
  gEnome = as.numeric(unname(gFf[gFf$type=="region", c("start", "end")]))
  gEnome = gEnome[order(gEnome)]

  ## transposon titre initiation
  tAg = "transposon population size per genome sd"
  tPn.pop = round(1/rNumVec(
    f = pMs$Value[pMs$Type==sub(" sd", " distribution", tAg)],
    L = pMs$Value[pMs$Type=="host organism constant population size"],
    p1 = pMs$Value[pMs$Type==sub(" sd", " mean", tAg)],
    p2 = pMs$Value[pMs$Type==tAg]),0)
  return(list(params = pMs, gene = gEne, transposon.titre = tPn.pop, genome = sum(c(gEne$length,gEne$interLength))))
}

##### Initiate genome pool #####
ini.host = function(host.var, gene.df, gene.var){
  nGene = nrow(gene.df)/2
  if(length(grep(";", gene.var))>0){gene.var = strsplit(gene.var, ";")[[1]]}
  gene.var = as.numeric(gene.var)
  if(any(gene.var > length(LETTERS))){
    gene.var[gene.var > length(LETTERS)] = length(LETTERS)
    warning(paste0("Too many alleles stated for genes ",paste(which(gene.var > length(LETTERS)), collapse = ","), " reset to ", length(LETTERS), ".\n"))
  }
  if(any(gene.var < 1)){
    gene.var[gene.var < 1] = 1
    warning(paste0("Too few alleles stated for genes ",paste(which(gene.var < 1), collapse = ","), " reset to 1.\n"))
  }
  if(!(length(gene.var) %in% c(1,nGene))){
    if(length(gene.var) > nGene){
      warning(paste0("Too many notations, extra ones are ignored. You provided ",length(gene.var), " gene variation notations, your gff file only has ",nGene, " genes."))
      gene.var = gene.var[1:nGene]
    }else{
      stop(paste0("Either one variation or specify variations for each gene in the genome. You only provided ",length(gene.var), " gene variation notations, your gff file has ",nGene, " genes."))
    }
  }else if(length(gene.var)==1){gene.var = rep(gene.var, nGene)}
  hOst = matrix(NA,nrow = as.numeric(host.var), ncol = nGene)
  for(i in 1:ncol(hOst)){ hOst[,i] = sample(LETTERS[1:gene.var[i]], size = host.var, replace = T)};rm(i)
  rEcessive = matrix(runif(prod(dim(hOst))), ncol = ncol(hOst))>.5
  hOst[rEcessive] = tolower(hOst[rEcessive])
  hOst.compact = apply(hOst, 1, function(x){paste0(x, collapse = "")})
  return(hOst.compact)
}

##### Rescale random number #####
reZero = function(x, new0 = 0, new1 = 1){return((x-new0)/(new1-new0))}

##### Transposon data format conversion #####
tPn.io = function(x, ref = tPn.0){
  if(class(x)=="data.frame"){
    x$paste = apply(x,1,function(x0){paste(x0, collapse = "!")})
    return(gsub(" ", "", paste(x$paste, collapse = ";")))
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
  if(length(x)==ncol(ref)){return(gsub(" ", "", paste(x, collapse = "!")))}
  stop("Provided ",length(x)," value(s) but ",ncol(ref)," values are needed: ",paste(colnames(ref), collapse = ", "),".")
}

##### Reconstruct gene table (one host) #####
reGeneDF = function(tPn, gene.df){
  if(tPn != ""){
    tPn = strsplit(tPn, ";")[[1]]
    for(i in 1:length(tPn)){
      i0 = tPn.io(tPn[i])
      i1 = which(gene.df$locus_tag == substr(i0[1],2,nchar(i0[1])))
      if(substr(i0[1],1,1) == "g"){
        gene.df$length[i1] = gene.df$length[i1] + as.numeric(i0[6])
        if(i1<nrow(gene.df)){
          gene.df$start[(i1+1):nrow(gene.df)] = gene.df$start[(i1+1):nrow(gene.df)] + as.numeric(i0[6])
        }
      }else{
        gene.df$interLength[i1] = gene.df$interLength[i1] + as.numeric(i0[6])
        gene.df$start[i1:nrow(gene.df)] = gene.df$start[i1:nrow(gene.df)] + as.numeric(i0[6])
      }
      gene.df$end[i1:nrow(gene.df)] = gene.df$end[i1:nrow(gene.df)] + as.numeric(i0[6])
  };rm(i)}
  return(gene.df)
}

##### New host population #####
host.reproduce = function(res.pool, gene.df, fitness.advantage, cell = "haploid"){
  # res.pool: 2 columns - $host, host genomes; $transposon, transposon notations
  ## Calculate ecological fitness deficit
  fitness.advantage = as.numeric(fitness.advantage)
  res.tmp = data.frame(host = unlist(read.table(text = res.pool$host, sep = ";")), transposon = NA, familyTree = unlist(read.table(text = res.pool$familyTree, sep = ";")), offspring.prob = 1)
  for(i in 1:nrow(res.pool)){if(length(grep(";",res.pool$transposon[i]))>0){
    t.tmp = tPn.io(res.pool$transposon[i])
    t.tmp$nonessential = as.numeric(!(t.tmp$gene %in% paste0("g",gene.df$locus_tag[gene.df$essential])))
    t.tmp$advantage = fitness.advantage/100 * as.numeric(t.tmp$gene %in% paste0("g",gene.df$locus_tag[gene.df$advantage]))
    t.tmp$sizeeff = t.tmp$size/sum(gene.df$length + gene.df$interLength)*-2
    g0 = substr(t.tmp$gene,2,2)==0

    if(sum(g0)>0){
      res.tmp$offspring.prob[i] = res.tmp$offspring.prob[i] * prod(t.tmp$nonessential[g0]) * (1 + sum(t.tmp$advantage[g0] + t.tmp$sizeeff[g0]))
      res.tmp$transposon[i] = tPn.io(t.tmp[g0,1:(ncol(t.tmp)-3)])
    }
    if(sum(g0)<length(g0)){
      res.tmp$offspring.prob[nrow(res.pool) + i] = res.tmp$offspring.prob[nrow(res.pool) + i] * prod(t.tmp$nonessential[!g0]) * (1 + sum(t.tmp$advantage[!g0] + t.tmp$sizeeff[!g0]))
      res.tmp$transposon[nrow(res.pool) + i] = tPn.io(t.tmp[!g0,1:(ncol(t.tmp)-3)])
    }
    rm(t.tmp,g0)
  }else if(length(grep("!",res.pool$transposon[i]))>0){
    t.tmp = tPn.io(res.pool$transposon[i])
    i0 = i + nrow(res.pool) * as.numeric(substr(t.tmp["gene"],2,2))
    res.tmp$offspring.prob[i0] = res.tmp$offspring.prob[i0] * !(t.tmp["gene"] %in% paste0("g",gene.df$locus_tag[gene.df$essential])) * (1 + sum((t.tmp["gene"] %in% paste0("g",gene.df$locus_tag[gene.df$advantage])) + as.numeric(t.tmp["size"])/sum(gene.df$length + gene.df$interLength)*-2))
    res.tmp$transposon[i0] = res.pool$transposon[i]
    rm(i0)
  }};rm(i)
  res.tmp$offspring.prob = res.tmp$offspring.prob/sum(res.tmp$offspring.prob)

  ## Sprouting offspring
  if(cell=="haploid"){
    offspring = sample(1:nrow(res.tmp), nrow(res.pool), replace = T, prob = res.tmp$offspring.prob)
    o.tmp = data.frame(host = res.tmp$host[offspring], transposon = gsub("i1","i0", gsub("g1","g0", res.tmp$transposon[offspring])), familyTree = ifelse((offspring %% nrow(res.pool))==0, nrow(res.pool), (offspring %% nrow(res.pool))))
    o.tmp[is.na(o.tmp)] = ""
    o.tmp = data.frame(host = paste(o.tmp$host,o.tmp$host, sep = ";"), transposon = paste(o.tmp$transposon,gsub("i0","i1", gsub("g0","g1", o.tmp$transposon)), sep = ";"), familyTree = paste(o.tmp$familyTree,o.tmp$familyTree, sep = ";"))
    if(length(grep("^;", o.tmp$transposon))>0){o.tmp$transposon[grep("^;", o.tmp$transposon)] = ""}
  }else{
    offspring = matrix(sample(1:nrow(res.tmp), nrow(res.pool)*2, replace = T, prob = res.tmp$offspring.prob), ncol = 2)
    o.tmp = as.data.frame(matrix(nrow = nrow(res.pool), ncol = ncol(res.pool)))
    colnames(o.tmp) = colnames(res.pool)
    o.tmp$host = paste(res.tmp$host[offspring[,1]],res.tmp$host[offspring[,2]], sep = ";")
    o.tmp$transposon = gsub(";$", "", gsub("^;", "", paste(gsub("i1","i0", gsub("g1","g0", res.tmp$transposon[offspring[,1]])), gsub("i0","i1", gsub("g0","g1", res.tmp$transposon[offspring[,2]])), sep = ";")))
    offspring = offspring %% nrow(res.pool)
    offspring[offspring==0] = nrow(res.pool)
    o.tmp$familyTree = paste(offspring[,1], offspring[,2], sep = ";")
  }
  return(o.tmp)
}

##### Modify notations according to transposon overlaps #####
tPn.x = function(tPn1, tPn2){
  t1 = tPn.io(tPn1); t2 = tPn.io(tPn2)
  if(t1[1] == t2[1]){
    t1.span = as.numeric(t1[2])+(1:as.numeric(t1[6]))-1
    if(as.logical(t1[4]) & (as.numeric(t1[2]) < as.numeric(t2[2])) & (as.numeric(t2[2]) %in% t1.span)){ # tPn1 assume always insert earlier than tPn2
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
  return(apply(tPn, 1, function(x){gsub(" ", "", paste0(x, collapse = "!"))}))
}

##### Single transposon relocation #####
tPn.reloc = function(tPn, gen, gene.df, tPn.prob=2){
  tPn = tPn.io(tPn)
  if(tPn.prob == 2){tPn.prob = runif(1)}
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
    h1.mod(pRobs[1], tPn.bg, sCene[1], ifelse(gToxic, ifelse(gTx[1], 0, ifelse(gTx[2], pAram[3], 1)),1)), # jump
    h1.mod(pRobs[2], tPn.bg, sCene[2], 1) # copy is not affected by genotoxic effect
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
tPn.act = function(tPn, equivalent, gen, gene.df, pAram, gToxic = F){ # gene.df must have all transposons attached
  x = tPn.io(tPn); x0 = tPn.io(tPn.io(tPn))
  x.probs = sCene.mod(pRobs = as.numeric(x[c("jumpRate", "copyRate")]), tPn.bg = tPn.get(gene.df), sCene = x[c("jumpH1", "copyH1")], pAram = pAram, gToxic = gToxic)
  x[c("jumpRate", "copyRate")] = x.probs[3:4]
  a = c(rNumVec(f = "uniform", L = 2, p1 = 0, p2 = 1) < x.probs[1:2], max(1, rpois(1, rnorm(1, 1, .01)))) # jump?, copy?, num copies
  a[1] = ifelse(gen > 0,a[1],1);a[2] = ifelse(gen > 0,a[2],0)
  colnames(gene.df)[1] = tPn.get(gene.df, F)
  if(as.logical(x["valid"])){
    if(a[1]>0){ # jump
      gene.df = tPn.reloc(tPn, gen, gene.df)
      x0 = tPn.get(gene.df)
      if(equivalent){x0 = c(tPn, x0)}
      colnames(gene.df)[1] = tPn.get(gene.df, F)
    }
    if(a[2]>0){ # copy
      x0 = unique(c(tPn,x0))
      for(i in 1:a[3]){
      gene.df = tPn.reloc(tPn, gen, gene.df)
      tTpn = c(tPn, tPn.get(gene.df))
      tCheck = abs(unlist(apply(tPn.io(paste(tTpn, collapse = ";")), 1, function(x){which(gene.df[,1]==substr(x[1], 2, nchar(x[1])))})) - nrow(gene.df)/2)
      if(x["copyDir"] == "terminus"){ # assume transposon is not quick enough to copy into the adjacent intragenic-gene pair
        cpValid = tCheck[2] < tCheck[1]
      }else if(x["copyDir"] == "origin"){
        cpValid = tCheck[2] > tCheck[1]
      }else{ cpValid = T }
      if(cpValid){
        x0 = c(x0, tPn.get(gene.df))
      }
      colnames(gene.df)[1] = tPn.get(gene.df, F)
      }}
  }else{x0 = tPn}
  colnames(gene.df)[1] = paste0(c(x0, colnames(gene.df)[1]), collapse = ";")
  return(gene.df)
}

##### Single gene recombination #####
gene.recom = function(h1.t, h2.t, g2to1, locusTags){
  if(h1.t=="" & h2.t==""){return(h1.t)}
  if(length(grep(paste0("g", locusTags[g2to1]), h1.t)) > 0){
    h1.t = strsplit(h1.t, ";")[[1]]
    h1.t = h1.t[-grep(paste0("g", locusTags[g2to1]), h1.t)]
    h1.t = paste0(h1.t, collapse = ";")
  }
  if(length(grep(paste0("g", locusTags[g2to1]), h2.t)) > 0){
    h2.t = strsplit(h2.t, ";")[[1]]
    h1.t = paste0(c(h1.t, h2.t[grep(paste0("g", locusTags[g2to1]), h2.t)]), collapse = ";")
  }
  return(h1.t)
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
    gene.df$recom = 0
    x.gene = sample(1:nrow(gene.df), numGenes[i], replace = F) # which gene being recombined (origins)
    gene.df$recom[x.gene] = 1
    gene.df$recomProb = runif(nrow(gene.df))
    for(i0 in x.gene){
      gene.df$recom = gene.df$recom + 2^(-abs(gene.df$start - gene.df$start[i0])/1000)
    };rm(i0)
    gene.df$recom = gene.df$recom > gene.df$recomProb

    ## remove intragenic transposons if consecutive genes are recombined (host independent)
    if(length(grep("!", res.pool$transposon[i]))>0){if((length(grep("i0", res.pool$transposon[i])) + length(grep("i1", res.pool$transposon[i])))>0){
      gene.df$consecutive = 0
      gene.df$consecutive[-nrow(gene.df)] = gene.df$recom[-nrow(gene.df)] + gene.df$recom[-1]
      g.bounds = c(range(which(substr(gene.df$locus_tag,1,1)==0)), range(which(substr(gene.df$locus_tag,1,1)==1)))
      gene.df$consecutive[g.bounds[1]] = gene.df$recom[g.bounds[1]] + gene.df$recom[g.bounds[2]]
      gene.df$consecutive[g.bounds[2]] = gene.df$recom[g.bounds[2]] + gene.df$recom[g.bounds[2]-1]
      gene.df$consecutive[g.bounds[3]] = gene.df$recom[g.bounds[3]] + gene.df$recom[g.bounds[4]]
      gene.df$consecutive[g.bounds[4]] = gene.df$recom[g.bounds[4]] + gene.df$recom[g.bounds[4]-1]
      t.tmp = tPn.io(res.pool$transposon[i])
      if(length(grep(";", res.pool$transposon[i]))>0){
        t.tmp$rm = t.tmp$gene %in% paste0("i",gene.df$locus_tag[gene.df$consecutive==2])
        res.pool$transposon[i] = tPn.io(t.tmp[!t.tmp$rm,-ncol(t.tmp)])
      }else{
        if(t.tmp["gene"] %in% paste0("i",gene.df$locus_tag[gene.df$consecutive==2])){res.pool$transposon[i] = ""}
    }}}

    ## Recombine genes
    x.recom = data.frame(
      gene = which(gene.df$recom),
      host = sample(c(1:nrow(res.pool))[-i], sum(gene.df$recom), replace = T)) # which host is the source
    x.recom$gene = x.recom$gene + ifelse(x.recom$gene > nrow(gene.df)/2, 1, 0)
    x.recom$gVar = substr(res.pool$host[x.recom$host],x.recom$gene,x.recom$gene)

    ## Reconstruct recipient genome
    h1.G = strsplit(res.pool$host[i], "")[[1]]
    h1.G[x.recom$gene] = x.recom$gVar
    res.pool$host[i] = paste0(h1.G, collapse = "")

    ## Reconstruct recipient transposons
    for(i0 in 1:nrow(x.recom)){
      res.pool$transposon[i] = sub(";$","",sub("^;","",gene.recom(h1.t = res.pool$transposon[i], h2.t = res.pool$transposon[x.recom$host[i0]], g2to1 = x.recom$gene[i0], locusTags = gene.df$locus_tag)))
    };rm(i0)
  }};rm(i)
  return(res.pool)
}
