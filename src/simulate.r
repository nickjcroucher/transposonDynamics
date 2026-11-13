#!/bin/env Rscript
# author: ph-u
# script: simulate.r
# desc: Discrete-time model for transposon ecology within one genome
# in: Rscript simulate.r [../custom/loc/input.csv] [../custom/loc/seed.csv] [seed number] [../custom/loc/scenario.csv] [scenario number]
# out: NA
# arg: 5 (optional)
# date: 20260623, 20261006

argv=(commandArgs(T))
if(length(argv) < 5){argv=c("../raw/input.csv", "../raw/seed.csv", "1", "../raw/scenario.csv", "186661")}
eLim = 1 # tag for transposon survival

##### env set #####
cat(date(),": set environment",argv[3],"-",argv[5],"\n")
source("func.r")
set.seed(read.csv(argv[2], header = F)[,1][as.numeric(argv[3])])
geneIn = gIn()
recomG = geneIn$gene$locus_tag[!is.na(geneIn$gene$recombination)]
rG = unique(substr(recomG, 2, nchar(recomG)))

sCene = read.csv(argv[4], header = FALSE, skip = as.numeric(argv[5]), nrows = 1, col.names = c("transposon", "host"), colClasses = c("character", "integer"), na.strings = "")
tPn.0 = tPn.0[which(tPn.0$uniqID %in% strsplit(sCene$transposon, ";")[[1]]),] # allow multiple transposons in simulation
host.0 = host.0[unique(sCene$host),]
gEn.max = xIn("host organism constant generation number")
perturbGen = ceiling(rev(seq(1, gEn.max, gEn.max/(host.0$genotoxic+1))))-1
perturbGen[perturbGen==0] = gEn.max + 1 # a numeric placeholder that can never achieve

##### Initiate genomic variation in populations #####
cat(date(),": initiate population",argv[3],"-",argv[5],"\n")
hOst = ini.host(geneIn$numGene, xIn("host genome variation"), xIn("host genetic variation")) # template for reconstruction
tPn = data.frame(ini = strsplit(tPn.io(tPn.0), ";")[[1]], uniqID = tPn.0$uniqID)

##### Initiate record dataframes #####
rec.host = rec.transposon = rec.offspring = rec.gRecom = rec.tpnCount = matrix(NA, nrow = gEn.max + 1, ncol = hPop <- xIn("host organism constant population size"))
rec.Fitness = cbind(rec.host,rec.host); rec.Fitness[1,] = 1/(hPop*2)
sim.df.col = c("familyTree", "hARecomRate", "transposon", "tpnFail", "host")
sim.df = as.data.frame(matrix("", nrow = hPop, ncol = length(sim.df.col)))
colnames(sim.df) = sim.df.col

##### Initial population #####
# Transposons (duplicating between first and second genomes)
ini.tPn = data.frame(tpn = sample(tPn$ini, sum(geneIn$transposon.titre), replace = T), loc = runif(n = sum(geneIn$transposon.titre), max = .5), host = rep(seq_len(length(geneIn$transposon.titre)), geneIn$transposon.titre), newLoc = NA)

## 1 Check no essential genes are hit
g0 = paste0("g", geneIn$gene$locus_tag[which(geneIn$gene$fitness == 0)])
repeat{
  ini.tPn = tPn.newloc(ini.tPn)
  ini.tPn$gene = tPn.io(ini.tPn$newLoc)$gene %in% g0
  if(all(!ini.tPn$gene)){break}
  ini.tPn$loc[ini.tPn$gene] = runif(n = sum(ini.tPn$gene), max = .5)
};rm(g0)

## 2 Check overlaps
ini.tPn = aggregate(newLoc ~ host, data = ini.tPn, FUN = function(x){paste0(x, collapse = ";")})
ini.tPn$newLoc = apply(ini.tPn, 1, FUN = function(x){return(tPn.x(x[2]))})

# Genotypes
fTree = as.data.frame(matrix(rep(sample(seq_len(length(hOst)), hPop * ifelse(host.0$cell=="haploid", 1, 2), replace = T), ifelse(host.0$cell=="haploid", 2, 1)), ncol = 2))
fTree$V1.tpn = ini.tPn$newLoc #[match(fTree$V1, ini.tPn$host)] # (transposon according to genotype)
fTree$V2.tpn = tPn.flat(x = ini.tPn$newLoc, ori = 0) #[match(fTree$V2, ini.tPn$host)]
sim.df$familyTree = paste(fTree[,1], fTree[,2], sep = ";")
sim.df$hARecomRate = host.0$homologousAutoRecom
sim.df$transposon = paste(fTree[,3], fTree[,4], sep = ";")
rm(ini.tPn)

##### Wright-Fisher / Neutral model run #####
gEn = 0; t.ori = Sys.time()
repeat{
  eTa(t.ori, gEn + 1, gEn.max + 1, paste("run simulation",argv[3],"-",argv[5]))
  # Population dynamics snapshot
  rec.host[gEn + 1,] = sim.df$host
  rec.transposon[gEn + 1,] = sim.df$transposon
  rec.offspring[gEn + 1,] = sim.df$familyTree
  rec.gRecom[gEn + 1,] = sim.df$hARecomRate
  rec.tpnCount[gEn + 1,] = lengths(strsplit(sim.df$transposon, ";"))
  save(hOst, rec.host, rec.transposon, rec.offspring, rec.gRecom, rec.tpnCount, rec.Fitness, sim.df, file = paste0("../data/tPn--", argv[3], "_", argv[5], ".rda"), compress = "xz")
  cat(" Analyzing;")
  gEn = gEn + 1; if(gEn > 1){source("sim_summary.r")}

  if(gEn > gEn.max){ break } # simulation done

  ##### Population reproduction stage #####
  cat(" Reproduction;")
  # 1 Reconstruct transposon-focused df
  s0.tpn = strsplit(gsub("^;$","", sim.df$transposon), ";")
  s0.tpn = data.frame(tpn = unlist(s0.tpn), host = rep(seq_len(length(s0.tpn)), lengths(s0.tpn)), fitness = 1)
  s0.tpn = cbind(s0.tpn[s0.tpn$tpn!="",], tPn.io(s0.tpn$tpn))
  s0.tpn$gi = substr(s0.tpn$tpn, 1, 1)
  s0.tpn$gene = substr(s0.tpn$gene, 2, nchar(s0.tpn$gene))
  s0.tpn$size = s0.tpn$size / geneIn$genome[2]
  s0.tpn$chromosome = as.numeric(substr(s0.tpn$gene,1,1))

  # 2 Calculate gene fitness
  s0.tpn$fitness = geneIn$gene$fitness[match(s0.tpn$gene, geneIn$gene$locus_tag)]
  s0.tpn$fitness[s0.tpn$gi=="i"] = 1 # assume intergenic regions do not affect fitness
  s0.tpn$fitness = s0.tpn$fitness - s0.tpn$size * host.0$transposonEffect
  s0.tpn$fitness = ifelse(s0.tpn$fitness < 0, 0, s0.tpn$fitness)

  # 3 Aggregate host fitness
  s0.fit = aggregate(fitness ~ host, data = s0.tpn, FUN = prod)
  if(nrow(s0.fit) != hPop){
    hPop.0 = seq_len(hPop)[!(seq_len(hPop) %in% s0.fit$host)]
    s0.fit.0 = as.data.frame(matrix(c(hPop.0, rep(1, length(hPop.0))), ncol = 2))
    colnames(s0.fit.0) = colnames(s0.fit)
    s0.fit = rbind(s0.fit, s0.fit.0)
    s0.fit = s0.fit[order(s0.fit$host),]
  }
  s0.fit$fitness = s0.fit$fitness * as.numeric(sim.df$tpnFail == "")
  if(sum(s0.fit$fitness)==0){ eLim = -2; break }

  # 4 Chromosome segregation
  s0.fit = rbind(s0.fit, s0.fit)
  s0.fit$chromosome = rep(c(0,1), each = hPop)
  s0.tpnSeg = aggregate(tpn ~ host + chromosome, data = s0.tpn, FUN = function(x){return(paste(x, collapse = ";"))})
  s0.tpnSeg = merge(s0.tpnSeg, s0.fit, all = T)
  s0.tpnSeg = s0.tpnSeg[order(s0.tpnSeg$chromosome, s0.tpnSeg$host),]
  s0.tpnSeg$tpn[is.na(s0.tpnSeg$tpn)] = ""
  s0.tpnSeg$tpn = tPn.flat(s0.tpnSeg$tpn, ori = 1)

  # 5 Repopulate
  rec.Fitness[gEn + 1,] = s0.tpnSeg$fitness/sum(s0.tpnSeg$fitness)
  fTree = as.data.frame(matrix(rep(sample(seq_len(nrow(s0.tpnSeg)), hPop * ifelse(host.0$cell=="haploid", 1, 2), replace = T, prob = as.numeric(rec.Fitness[gEn + 1,])), ifelse(host.0$cell=="haploid", 2, 1)), ncol = 2))
  fTree$V1.tpn = s0.tpnSeg$tpn[fTree$V1]
  fTree$V2.tpn = tPn.flat(s0.tpnSeg$tpn[fTree$V2], ori = 0)
  fTree$har = h1.mod((sim.df$hARecomRate[s0.tpnSeg$host[fTree$V1]] + sim.df$hARecomRate[s0.tpnSeg$host[fTree$V2]]) / 2, host.0$hAReomH1)
  # 6 Inherit appropriate parameters
  sim.df$familyTree = paste(fTree[,1], fTree[,2], sep = ";")
  sim.df$transposon = paste(fTree[,3], fTree[,4], sep = ";")
  sim.df$transposon[sim.df$transposon == ";"] = ""
  if(host.0$recom == 0 & any(sim.df$transposon=="")){stop("Transposon error on reproduction step")}
  sim.df$hARecomRate = fTree$har
  sim.df[,4:5] = ""

  ##### Transposon activating stage #####
  cat(" Transposons;")
  # 1 Set transposon players
  s0.tpn = strsplit(gsub("^;$","", sim.df$transposon), ";")
  s0.tpnLen = lengths(s0.tpn)
  if(sum(s0.tpnLen)==0){ eLim = -1; break }
  s0.tpn = data.frame(tpn = unlist(s0.tpn), host = rep(seq_len(length(s0.tpn)), s0.tpnLen), fitness = 1)
  if(host.0$recom == 0 & any(!(seq_len(hPop) %in% s0.tpn$host))){stop("Transposon error on transposon expansion step")}
  s0.tpn = cbind(s0.tpn[s0.tpn$tpn!="",], tPn.io(s0.tpn$tpn))
  s0.tpn$numTpn = rep(s0.tpnLen, s0.tpnLen)

  # 2 Activity rate modifications
  s0.tpn$jumpRateFinal = s0.tpn$jumpRate
  s0.tpn$copyRateFinal = s0.tpn$copyRate

  ## 2.1 H1: Charlesworth and/or Evolving?
  s0.tpn$jumpRateFinal = as.numeric(apply(s0.tpn, 1, FUN = function(x){h1.mod(ori = x[16], h1 = x[11], numTpn = x[15])}))
  if(length(grep("evolv", s0.tpn$jumpH1)) > 0){s0.tpn$jumpRate[grep("evolv", s0.tpn$jumpH1)] = s0.tpn$jumpRateFinal[grep("evolv", s0.tpn$jumpH1)]}

  s0.tpn$copyRateFinal = as.numeric(apply(s0.tpn, 1, FUN = function(x){h1.mod(ori = x[17], h1 = x[13], numTpn = x[15])}))
  if(length(grep("evolv", s0.tpn$copyH1)) > 0){s0.tpn$copyRate[grep("evolv", s0.tpn$copyH1)] = s0.tpn$copyRateFinal[grep("evolv", s0.tpn$copyH1)]}

  ## 2.2 genotoxic?
  if(gEn %in% perturbGen){
    s0.tpn$jumpRateFinal = s0.tpn$jumpRateFinal * ifelse(runif(nrow(s0.tpn)) > xIn("percentage transposon perturbation genotoxic")/100, 1, 0) * ifelse(runif(nrow(s0.tpn)) > xIn("percentage chance transposon perturbation boost")/100, 1, xIn("percentage amplitude transposon perturbation boost")/100)
  }

  # 3 Valid for activity? Jump? Copy? Replicating number? Replacement if jumped?
  ## 3.1 Set probabilities for activities
  s0.tpn[,c("mv", "mvhAR", "mvLoc", "cp", "cpRep")] = matrix(c(runif(nrow(s0.tpn) * 4), rpois(nrow(s0.tpn), 1)), ncol = 5)
  while(any(s0.tpn$cpRep < 1)){s0.tpn$cpRep[which(s0.tpn$cpRep < 1)] = rpois(sum(s0.tpn$cpRep < 1), 1)}
  s0.tpn$mv = ifelse(s0.tpn$valid, ifelse(s0.tpn$mv > s0.tpn$jumpRateFinal, NA, s0.tpn$mv), NA)
  s0.tpn$mvhAR = !(s0.tpn$mvhAR > as.numeric(sim.df$hARecomRate[s0.tpn$host])) # TRUE = success homologous recombine
  s0.tpn$cp = ifelse(s0.tpn$valid, ifelse(s0.tpn$cp > s0.tpn$copyRateFinal, NA, s0.tpn$cp), NA)

  ## 3.2 Homolog lookup
  s0.tpn$mvSub = NA
  for(i in unique(s0.tpn$host)){
    r = which(s0.tpn$host == i)
    hom = paste0(substr(s0.tpn$gene[r],1,1), (as.numeric(substr(s0.tpn$gene[r],2,2)) + 1) %% 2, substring(s0.tpn$gene[r],3))
    s0.tpn$mvSub[r] = r[match(paste(hom, s0.tpn$location[r]), paste(s0.tpn$gene[r], s0.tpn$location[r]))]
  };rm(i, r, hom)

  # 4 Scenario deductions
  ## 4.1 tpnFail
  tFail = s0.tpn[which(!is.na(s0.tpn$mv) & !s0.tpn$mvhAR), c("gene", "host")]
  for(i in unique(tFail$host)){sim.df$tpnFail[i] = paste0(tFail$gene[which(tFail$host==i)], collapse = ";")};rm(i, tFail)

  ## 4.2 Potentially-modified parameter incorporation
  s0.tpn$tpnEvo = strsplit(tPn.io(s0.tpn[,which(colnames(s0.tpn)==colnames(tPn.0)[1]):which(colnames(s0.tpn)==colnames(tPn.0)[ncol(tPn.0)])]), ";")[[1]]

  ## 4.3 Identify lethal double strand breaks (prokaryotes only)
  if(host.0$cell == "haploid"){
    tFail = s0.tpn[which(!is.na(s0.tpn$mv)),]
    tFail = unique(tFail$host[which(row.names(tFail) %in% tFail$mvSub)])
    sim.df$tpnFail[tFail] = sub("^;","",paste(sim.df$tpnFail[tFail],"@",sep = ";"))
  }

  # 5 Transposons migration
  i.mv = which(!is.na(s0.tpn$mv))
  i.cp = which(!is.na(s0.tpn$cp))

  ## 5.1 mv
  if(length(i.mv) > 0){
    s0.final = s0.tpn[-i.mv, 1:3]
    s0.sub = s0.tpn[as.numeric(s0.tpn$mvSub[i.mv[s0.tpn$mvhAR[i.mv] & !is.na(s0.tpn$mvSub[i.mv])]]),1:3] # jump tpn replacements
    s0.sub$tpn = paste0(substr(s0.sub$tpn,1,1), (as.numeric(substr(s0.sub$tpn,2,2)) + 1) %% 2, substr(s0.sub$tpn,3,nchar(s0.sub$tpn))) # chromosome location change
    tpns = s0.tpn[i.mv,1:3]
    tpns$tpn = s0.tpn$tpnEvo[i.mv]
    tpns$newLoc = runif(nrow(tpns)) # draw new locations
    for(i in unique(tpns$host)){ # record new locations into tpn notations
      i0 = reGeneDF(tPn = paste0(s0.tpn$tpn[s0.tpn$host==i], collapse = ";"), gene.df = geneIn$gene)
      for(i1 in which(tpns$host==i)){
        i2 = tPn.reloc(num = tpns$newLoc[i1], tpn = tpns$tpn[i1], gEn = gEn, gDF = i0)
        tpns$tpn[i1] = i2$tpn
        i0 = i2$gDF
      };rm(i1)
    };rm(i,i0,i.mv)
    s0.final = rbind(s0.final, tpns[,1:3])
    s0.final = rbind(s0.final, s0.sub); rm(s0.sub)
  }else{ s0.final = s0.tpn[, 1:3] }

  ## 5.2 cp
  if(length(i.cp) > 0){
    tpns = s0.tpn[rep(i.cp, s0.tpn$cpRep[i.cp]),]
    tpns$tpn = tpns$tpnEvo
    tpns$mvLoc = runif(nrow(tpns))
    for(i in unique(tpns$host)){ # record new locations into tpn notations
      i0 = reGeneDF(tPn = paste0(s0.tpn$tpn[s0.tpn$host==i], collapse = ";"), gene.df = geneIn$gene)
      for(i1 in which(tpns$host==i)){
        i2 = tPn.reloc(num = tpns$mvLoc[i1], tpn = tpns$tpn[i1], gEn = gEn, gDF = i0)
        if(tpns$copyDir[i1] == "both"){i3 = c(tpns$copyDir[i1], "terminus", "origin")}else{i3 = tpns$copyDir[i1]}
        if(tPn.comLoc(oldNew = c(tpns$tpn[i1], i2$tpn), gDF = i0) %in% i3){
          tpns$tpn[i1] = i2$tpn
          i0 = i2$gDF
      }else{ tpns$tpn[i1] = NA }};rm(i1, i3)
    };rm(i,i0,i.cp)
    s0.final = rbind(s0.final, tpns[!is.na(tpns$tpn),1:3])
  }

  ## 5.3 validity check
  for(i in unique(s0.final$host)){
    r = which(s0.final$host == i)
    s0.final$tpn[r] = strsplit(tPn.x(paste0(s0.final$tpn[r], collapse = ";")), ";")[[1]]
  };rm(i, r)
  if(nrow(s0.final)<1){ eLim = -1; break }
  s0.final = cbind(s0.final, tPn.io(s0.final$tpn))
  if(host.0$recom == 0 & (length(seq_len(hPop)) != length(unique(s0.final$host)))){stop("Transposon error at end of transposon stage")}

  ##### Gene recombination stage (assume no transposons excise / add complications) #####
  cat(" Recombination\n")
  # 1 Check if there is recombination happening
  reCom = host.0$recom > 0
  if(reCom){
    reCom.df = as.data.frame(matrix(0, nrow = hPop, ncol = length(rG)))
    colnames(reCom.df) = rG
    for(i in seq_len(length(rG))){
      i0 = s0.final[grep(paste0(rG[i], "!"), s0.final$tpn, fixed = TRUE), c("host", "gene")]
      i0 = table(unique(data.frame(host = i0$host, copy = substr(i0$gene, 2, 2)))$host) # disruption occur when transposon inserted onto the gene / intergenic region in front; copies disrupted per host
      if(length(i0) > 0){ reCom.df[as.numeric(names(i0)),i] = i0 }
    };rm(i,i0)

  # 2 Number of genes going to be recombined
    gRecRate = host.0$recom * nrow(geneIn$gene) / 5000 # assume recombination mean length = 5kb [202610231012]

  ## Genome size increase due to transposons
    numGenes = aggregate(size ~ host, data = s0.final, FUN = sum)
    if(nrow(numGenes) < hPop){
      i = which(!(seq_len(hPop) %in% s0.final$host))
      i = as.data.frame(matrix(c(i, rep(0, length(i))), ncol = 2))
      colnames(i) = colnames(numGenes)
      numGenes = rbind(numGenes, i); rm(i)
      numGenes = numGenes[order(numGenes[,1]),]
    }

  ## 2.2 Number of recombined genes
    numGenes = gRecRate * (mean(geneIn$gene$length) * (1 + numGenes[,2] / geneIn$genome[2] * (host.0$recomH1!="switch")) * as.numeric(!(apply(reCom.df, 1, max) == 2)))
    for(i in seq_len(hPop)){numGenes[i] = rpois(1, lambda = numGenes[i])}
    numGenes = pmin(numGenes, nrow(geneIn$gene))
    if(all(apply(reCom.df, 1, max) == 2) | all(numGenes==0)){reCom = FALSE}
  }

  if(reCom){
  # 3 Hotspots selection
  rGenes = c(); for(i in seq_len(hPop)){rGenes = c(rGenes, sample(geneIn$gene$locus_tag, numGenes[i], replace = FALSE))};rm(i)
  rGenes = data.frame(gene = rGenes, host = rep(seq_len(hPop), numGenes))
  hRecom = which(numGenes>0)

  # 4 Nearby gene regions
  rGSel = NULL; for(i in hRecom){ rGSel = rbind(rGSel, cbind(gRc.prob(rGenes$gene[rGenes$host==i]), i)) };rm(i)
  colnames(rGSel)[3] = "host"

  # 5 Calculate the probability of each gene in each host being selected
  rGenes = matrix(rep(geneIn$gene$length, hPop), ncol = hPop)
  tpns = aggregate(size ~ host + gene, data = s0.final, FUN = sum)
  tpns = tpns[substr(tpns$gene,1,1)=="g",]
  tpns$gene = match(substr(tpns$gene,2,nchar(tpns$gene)), geneIn$gene$locus_tag)
  for(i in seq_len(nrow(tpns))){rGenes[tpns$gene[i], tpns$host[i]] = rGenes[tpns$gene[i], tpns$host[i]] + tpns$size[i]};rm(i)
  rGenes = rGenes / rowSums(rGenes)

  ## Exclude broken genes in tpnFail being selected to recombine
  if(any(sim.df$tpnFail!="")){
    s0.tpn = strsplit(sim.df$tpnFail, ";")
    for(i in which(lengths(s0.tpn)>0)){
      s0.tpn[[i]] = s0.tpn[[i]][substr(s0.tpn[[i]],1,1)=="g"]
      if(length(s0.tpn[[i]])>0){ rGenes[geneIn$gene$locus_tag %in% substr(s0.tpn[[i]],2,nchar(s0.tpn[[i]])),i] = 0 }
    };rm(i)
  }

  # 6 Gene recombination
  rGSel$paste = paste(rGSel$contiguous, rGSel$host, sep = ";")

  ## 6.1 Get genes for recombination
  if(host.0$cell=="diploid"){ # gene recombine with its own
    rGSel$src = rGSel$host; rGSel$chromosome = (as.numeric(substr(rGSel$locus_tag, 1, 1)) + 1) %% 2
  }else{ # gene recombine with other host cells
    x0 = rle(rGSel$paste)
    rGSel$src = rep(sample(seq_len(hPop), length(x0$lengths), replace = TRUE), x0$lengths)
    while(any(rGSel$host==rGSel$src)){ ## Get host cell excluding oneself
      x1 = which(rGSel$host==rGSel$src)
      x2 = rle(rGSel$paste[x1])
      rGSel$src[x1] = rep(sample(seq_len(hPop), length(x2$lengths), replace = TRUE), x2$lengths)
      rm(x1,x2)
    }
    rGSel$chromosome = rep(sample(0:1, length(x0$lengths), replace = TRUE), x0$lengths)
    rm(x0)
  }
  rGSel$donor = paste0(rGSel$chromosome, substring(rGSel$locus_tag, 2))

  ## 6.2 Get transposons in genes / intergenic regions recombine
  rGSel$gLoc = match(rGSel$locus_tag, geneIn$gene$locus_tag)
  rGSel$same = c(FALSE, rGSel$src[-nrow(rGSel)] == rGSel$src[-1])
  rGSel$ghost = paste0(ifelse(rGSel$same, paste0("i",rGSel$locus_tag,";"), ""),paste0("g",rGSel$locus_tag))
  rGSel$shost = paste0(ifelse(rGSel$same, paste0("i",rGSel$donor,";"), ""),paste0("g",rGSel$donor))

  ## Double-strand break rescue
  x0 = gsub(";$", "", gsub("^;", "", gsub("@", "", sim.df$tpnFail)))
  if(any(x0!="")){
    x1 = strsplit(sim.df$tpnFail, ";")
    for(i in which(x0!="")){
      sim.df$tpnFail[i] = paste0(x1[[i]][!(x1[[i]] %in% strsplit(paste0(rGSel$ghost[rGSel$host==i], collapse = ";"), ";")[[1]])], collapse = ";")
    };rm(i,x1)
  };rm(x0)

  ## 6.3 Transposons elmination
  i0 = paste(s0.final$host, s0.final$gene, sep = "!")
  i1 = paste(rGSel$host, rGSel$ghost, sep = "!")
  i2 = regexpr("!", i1)
  for(i in grep(";", i1)){
    i1[i] = gsub(";", paste0(";", substr(i1[i],1,i2[i])), i1[i])
  };rm(i)
  i1 = unlist(strsplit(i1, ";"))
  sDF.rm = which(i0 %in% i1)
  rm(i1,i2)

  ## 6.4 Transposons immigration
  i1 = paste(rGSel$host, substr(rGSel$locus_tag,1,1), rGSel$src, rGSel$shost, sep = "!")
  i2 = gregexpr("!", i1)
  for(i in grep(";", i1)){
    i1[i] = gsub(";", paste0(";", substr(i1[i],1,i2[[i]][3])), i1[i])
  };rm(i)
  i1 = read.table(text = unlist(strsplit(i1, ";")), sep = "!")
  i1$paste = paste(i1[,3], i1[,4], sep = "!")
  sDF.rp = s0.final[0,1:3]
  for(i in unique(i1[,1])){
    i3 = which(i0 %in% i1[which(i1[,1] == i),5])
    if(length(i3) > 0){
      i3 = s0.final[i3,]
      for(i4 in seq_len(nrow(i3))){
        i3$tpn[i4] = paste0(substr(i3$tpn[i4],1,1), unique(i1[which(i1[,1]==i & i1[,3]==i3$host[i4] & i1[,4]==i3$gene[i4]),2]), substr(i3$tpn[i4],3,nchar(i3$tpn[i4])))
      };rm(i4)
      i3$host = i
      sDF.rp = rbind(sDF.rp, i3[,1:3])
      }
  };rm(i,i3)
  rm(i0,i1,i2)
  s0.final = rbind(s0.final[!(seq_len(nrow(s0.final)) %in% sDF.rm),1:3], sDF.rp); rm(sDF.rm, sDF.rp)
  if(nrow(s0.final)<1){ eLim = -1; break }

  # 7 Record gene recombinations (rec.host)
  rGSel$paste = paste(rGSel$gLoc, rGSel$src, rGSel$chromosome, sep = "!")
  rGSel = aggregate(paste ~ host, data = rGSel, FUN = function(x){paste(x, collapse = ";")})
  sim.df$host[hRecom] = rGSel$paste
  }

  ##### Harmonize all transposon records from all stages #####
  s0.final = aggregate(tpn ~ host, data = s0.final, FUN = function(x){paste0(x, collapse = ";")})
  sim.df$transposon = s0.final[match(as.numeric(row.names(sim.df)), s0.final$host),2]
  if(host.0$recom == 0 & any(is.na(sim.df$transposon))){stop("Transposon error at recombination stage")}
  if(any(is.na(sim.df$transposon))){sim.df$transposon[is.na(sim.df$transposon)] = ""}
  sim.df$transposon = apply(sim.df, 1, FUN = function(x){return(tPn.x(x[3]))})
}

cat(date(),": printing warnings",argv[3],"-",argv[5],"\n")
print(warnings())

##### Simulation record export #####
if(!(tolower(xIn("full result")) %in% c("yes", "true"))){
  cat(date(),": simplified result export",argv[3],"-",argv[5],"\n")
  oUt = unique(ceiling(seq(1,nrow(rec.host), (nrow(rec.host)-1)/10)))
  rec.host = rec.host[oUt,]; rec.transposon = rec.transposon[oUt,]; rec.offspring = rec.offspring[oUt,]; rec.gRecom = rec.gRecom[oUt,]; rec.tpnCount = rec.tpnCount[oUt,]
  save(hOst, rec.host, rec.transposon, rec.offspring, rec.gRecom, rec.tpnCount, sim.df, file = paste0("../data/tPn--", argv[3], "_", argv[5], ".rda"), compress = "xz")
}

if(eLim == -1){tGt = "Transposons"}else if(eLim == -2){tGt = "Host population"}
if(eLim < 0){ tXt = paste0(tGt," eliminated: gen ",gEn,";") }else{tXt = ""}
cat(date(),":",tXt,"simulation completed",argv[3],"-",argv[5],"\n")
# c0 = rep(0,nrow(rec.transposon));for(i in seq_len(length(c0))){c0[i] = length(unlist(strsplit(as.character(rec.transposon[i,]), ";")))};rm(i);plot(x = seq_len(length(c0)), y = c0, cex = .1)
