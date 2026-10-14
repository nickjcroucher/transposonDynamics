#!/bin/env Rscript
# author: ph-u
# script: func.r
# desc: self-defined functions for transposon ecology simulation
# in: source("func.r")
# out: NA
# arg: NA
# date: 20260623, 20260928

# library(stringr)

##### Constants #####
rAw = c(dirname(argv[1]), basename(argv[1]))
f.all = c(list.files(rAw[1], full.names = TRUE), list.files("../data", pattern = "gff$", full.names = TRUE))
tPn.0 = read.csv(f.all[grepl("-tpn\\.csv$", f.all)], header = T, na.strings = "")
host.0 = read.csv(f.all[grepl("-host\\.csv$", f.all)], header = T)
pMs = read.csv(argv[1], header = T)

##### Extract input parameters #####
xIn = function(x, ref = pMs){
  z0 = ref$Value[ref$Type==x]
  if(suppressWarnings(!is.na(as.numeric(z0)))){z0 = as.numeric(z0)}
  return(z0)
}

##### Reorganize GFF file #####
gffClean = function(fLst = f.all){
  gffs = grep("gff$", fLst, value = TRUE); hit = gffs[grepl(xIn("ref genome"), basename(gffs), fixed = TRUE)]
  if(length(hit) != 1){stop("expected one GFF named after ", xIn("ref genome"), "; found ", length(hit))}
  gFf = read.table(hit, sep = "\t", header = F, quote = "")
  colnames(gFf) = c("seqname","source","type","start","end","score","strand","phase","attributes")
  x = strsplit(gFf$attributes, ";")
  x0 = unlist(strsplit(sub("=","=@",unlist(x)), "="))
  cNam = unique(x0[grep("@", x0, invert = T)])
  aTt = as.data.frame(matrix(nrow = length(x), ncol = length(cNam)))
  colnames(aTt) = cNam
  for(i0 in seq_len(length(x))){for(i1 in seq_len(length(x[[i0]]))){
    i2 = strsplit(x[[i0]][i1], "=")[[1]]
    aTt[i0,colnames(aTt)==i2[1]] = i2[2]
  }};rm(i0,i1)
  gFf = cbind(gFf[,-ncol(gFf)],aTt)
  return(gFf)
}

##### Gene details input #####
gIn = function(){
  gFf = unique(gffClean()[,c("locus_tag","type","start","end")])
  gEnome = as.numeric(gFf[which(gFf$type=="region")[1],-(1:2)])
  gFf = unique(gFf[which(gFf$type!="region" & !is.na(gFf$locus_tag)),-2])

  ## Gene lengths: [inter1]->[gene1]->[inter2]->[gene2]...
  gFf$length = gFf$end - (gFf$start - 1)
  gFf$interLength = c(gEnome[2] - gFf$end[nG <- nrow(gFf)] + gFf$start[1] - gEnome[1], gFf$start[-1] - (gFf$end[-nrow(gFf)] + 1))
  g0 = which(gFf$interLength < 0)
  gFf$length[g0-1] = gFf$length[g0-1] + gFf$interLength[g0]; gFf$interLength[g0] = 0

  ## Gene fitness
  g0 = read.csv(xIn("gene fitness"), header = TRUE)
  gFf$fitness = g0$fitness[match(gFf$locus_tag, g0$locus_tag)]
  gFf$fitness[is.na(gFf$fitness)] = 0 # assume essential genes cannot generate transposon library mutants

  ## Gene recombination machinery
  gFf$recombination = gFf$locus_tag %in% strsplit(xIn("genes for recombination mechanism"), ";")[[1]]

  ## Diploid chromosome in cell
  gFf = rbind(gFf, gFf)
  gFf$locus_tag = paste0(rep(c(0,1), each = nG), gFf$locus_tag)
  row.names(gFf) = NULL
  gFf$recombination = ifelse(gFf$recombination>0, ifelse(as.numeric(row.names(gFf))>nG, 1, 0), NA)

  ## Initial transposon copy number
  tPn.pop = rpois(n = xIn("host organism constant population size"), lambda = xIn("transposon population size per genome mean"))
  while(any(tPn.pop < 1)){ tPn.pop[tPn.pop < 1] = rpois(n = sum(tPn.pop < 1), lambda = xIn("transposon population size per genome mean")) }

  return(list(gene = gFf, transposon.titre = tPn.pop, numGene = nG, genome = gEnome))
}

##### Initiate genome pool #####
ini.host = function(nGenes, host.var, gene.var){
  gv = rep_len(pmax(pmin(as.numeric(strsplit(as.character(gene.var), ";")[[1]]), length(letters)), 1), nGenes)
  h0 = matrix(vapply(gv, function(v){sample(c(LETTERS[seq_len(v)], letters[seq_len(v)]), host.var, replace = T)}, character(host.var)), nrow = host.var)
  return(apply(h0, 1, paste0, collapse = ""))
}

##### Transposon data format conversion #####
tPn.io = function(x, ref = tPn.0){
  if(class(x)=="data.frame"){
    x$paste = apply(x,1,function(x0){paste(x0, collapse = "!")})
    return(gsub(" ", "", paste(x$paste, collapse = ";")))
  }else{
    return(read.table(text = gsub(";", "\n", paste0(x, collapse = ";"), fixed = T), sep = "!", header = F, quote = "", comment.char = "", col.names = colnames(ref)))
  }
}

##### Reconstruct gene table (one host, multiple transposons) #####
reGeneDF = function(tPn, gene.df){
  if(tPn==""){return(gene.df)}
  tPn = tPn.io(tPn)
  tPn$cOl = match(ifelse(substr(tPn$gene, 1, 1)=="g","length","interLength"), colnames(gene.df))
  tPn$rOw = match(substr(tPn$gene, 2, nchar(tPn$gene)), gene.df$locus_tag)
  t0 = aggregate(size ~ rOw + cOl, data = tPn, FUN = sum)
  t1 = unique(t0$cOl)
  for(i in seq_len(length(t1))){
    gene.df[t0$rOw[which(t0$cOl==t1[i])], t1[i]] = gene.df[t0$rOw[which(t0$cOl==t1[i])], t1[i]] + t0$size[which(t0$cOl==t1[i])]
  };rm(i)
  return(gene.df)
}

##### Relocating transposons #####
tPn.reloc = function(num, tpn, gEn = 0, gDF = geneIn$gene){
  tpn = tPn.io(tpn)
  tpn$generation = gEn
  gDF$totLn = cumsum(gDF$length + gDF$interLength)
  num = min(ceiling(num * gDF$totLn[nrow(gDF)]), gDF$totLn[nrow(gDF)])
  lOc = which(gDF$totLn >= num)[1]
  if(length(lOc)==0){lOc = nrow(gDF)}
  tpn$gene = gDF$locus_tag[lOc]
  num = gDF$length[lOc] + gDF$interLength[lOc] - gDF$totLn[lOc] + num
  gI = num > gDF$interLength[lOc]
  tpn$gene = paste0(ifelse(gI, "g", "i"), tpn$gene)
  tpn$location = num - ifelse(gI, gDF$interLength[lOc], 0)
  gDF[lOc,ifelse(gI, "length", "interLength")] = gDF[lOc,ifelse(gI, "length", "interLength")] + tpn$size
  gDF$totLn = NULL
  return(list(tpn = tPn.io(tpn), gDF = gDF))
}

##### Check transposons overlap #####
tPn.x = function(tpns){
  if(tpns==""){return("")}
  tpns = tPn.io(tpns)
  tpns$valid = TRUE
  if(any(table(tpns$gene)>1)){
    i0 = which(duplicated(tpns$gene))
    for(i in seq_len(length(i0))){
      i.prior = which(tpns$gene[seq_len(i0[i]-1)] == tpns$gene[i0[i]])
      for(i1 in seq_len(length(i.prior))){
        if(tpns$valid[i.prior[i1]]){
          tpns$valid[i.prior[i1]] = !(tpns$location[i0[i]] %in% (seq_len(tpns$size[i.prior[i1]]) - 1 + tpns$location[i.prior[i1]]))
        }
      }
    }
  }
  return(tPn.io(tpns))
}

##### Flatten transposon chromosomal locating notation #####
tPn.flat = function(x, ori = 0){
  ori = ori %% 2
  for(i in c("g", "i")){
    x = gsub(paste0(i, ori), paste0(i, (ori+1) %% 2), x)
  }
  return(x)
}

##### Rate modifier #####
h1.mod = function(ori, h1, numTpn = 1){
  ori = as.numeric(ori); numTpn = as.numeric(numTpn); oriLen = length(ori)
  if(length(numTpn)!=1 & length(numTpn)!=length(ori)){stop("Length of numTpn can either be 1 or length of ori")}
  ## evolving based on the suppressed activity background
  if(length(grep("charlesworth", h1)) > 0){ ori = ori / numTpn }
  if(length(grep("evolv", h1)) > 0){ ori = ori * abs(rnorm(oriLen, mean = 1, sd = .01))}
  return(ori)
}

##### Set transposon location #####
tPn.newloc = function(tpnDF, gen = 0, gDF = geneIn$gene){
  for(i in seq_len(nrow(tpnDF))){
    if(i==1 || tpnDF$host[i] != tpnDF$host[i-1]){gD = gDF}else{gD = g.tmp$gDF}
    g.tmp = tPn.reloc(num = tpnDF$loc[i], tpn = tpnDF$tpn[i], gEn = gen, gDF = gD)
    tpnDF$newLoc[i] = g.tmp$tpn
  }
  return(tpnDF)
}

##### Location comparison between new & old transposon #####
tPn.comLoc = function(oldNew, G = geneIn$genome[2], gDF = geneIn$gene){ # input two vectors
  if(length(oldNew)!=2){stop("Length-two vectors indicating the old and new transposons are needed")}
  t0 = tPn.io(paste(oldNew, collapse = ";"))
  if(t0$gene[1]==t0$gene[2] && t0$location[1]==t0$location[2]){return(unique(t0$copyDir))}
  r = match(substr(t0$gene, 2, nchar(t0$gene)), gDF$locus_tag)
  pos = gDF$start[r] - 1 + t0$location - ifelse(substr(t0$gene, 1, 1)=="i", gDF$interLength[r], 0)
  pos = (pos - 1) %% G + 1 # including intergenic segment between origin & gene 1
  d = pmin(pos, G - pos) # distance from origin (position 1)
  if(d[2] > d[1]){return("terminus")}else{return("origin")}
} # if jump back onto the exact same location, allow survival

##### Gene hotspots additive probability #####
gRc.prob = function(rGs, g.hotspot = geneIn$gene){
  g.hotspot$chromosome = as.numeric(substr(g.hotspot$locus_tag,1,1))
  g.hotspot$hotspot = runif(nrow(g.hotspot))
  g.hotspot$hotspot[g.hotspot$locus_tag %in% rGs] = 0
  g.hotspot$threshold = 0
  g.0 = g.hotspot[g.hotspot$chromosome==0,]; g.1 = g.hotspot[g.hotspot$chromosome==1,]
  r.0 = rGs[substr(rGs,1,1)==0]; r.1 = rGs[substr(rGs,1,1)==1]
  for(z in seq_len(length(r.0))){
    g.0$threshold = g.0$threshold + 2^(-abs(g.0$start - g.0$start[which(g.0$locus_tag==r.0[z])])/1000)
  };rm(z)
  for(z in seq_len(length(r.1))){
    g.1$threshold = g.1$threshold + 2^(-abs(g.1$start - g.1$start[which(g.1$locus_tag==r.1[z])])/1000)
  };rm(z)
  g.hotspot = rbind(g.0, g.1)
  return(c(g.0$locus_tag[g.0$hotspot <= g.0$threshold], g.1$locus_tag[g.1$hotspot <= g.1$threshold]))
}
