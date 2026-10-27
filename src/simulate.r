#!/bin/env Rscript
# author: ph-u
# script: simulate.r
# desc: Discrete-time model for transposon ecology within one genome
# in: Rscript simulate.r [../custom/loc/input.csv] [../custom/loc/seed.csv] [seed number] [../custom/loc/scenario.csv] [scenario number]
# out: NA
# arg: 5 (optional)
# date: 20260623, 20261006

argv=(commandArgs(T))
if(length(argv) < 5){argv=c("../data/input.csv", "../raw/seed.csv", "1", "../raw/scenario.csv", "211681")} # 208225 vs 211681
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
sim.df.col = c("familyTree", "hARecomRate", "transposon", "tpnFail", "host")
sim.df = as.data.frame(matrix("", nrow = hPop, ncol = length(sim.df.col)))
colnames(sim.df) = sim.df.col

##### Initial population #####
# Genotypes
fTree = matrix(rep(sample(seq_len(length(hOst)), hPop * ifelse(host.0$cell=="haploid", 1, 2), replace = T), ifelse(host.0$cell=="haploid", 2, 1)), ncol = 2)
sim.df$familyTree = paste(fTree[,1], fTree[,2], sep = ";")
sim.df$hARecomRate = host.0$homologousAutoRecom

# Transposons (duplicating between first and second genomes)
ini.tPn = data.frame(tpn = sample(tPn$ini, sum(geneIn$transposon.titre), replace = T), loc = runif(n = sum(geneIn$transposon.titre), max = .5), host = rep(seq_len(nrow(sim.df)), geneIn$transposon.titre), newLoc = NA)

# 1 set new locations
ini.tPn = tPn.newloc(ini.tPn)

# 2 Check overlaps
sim.df$transposon = aggregate(newLoc ~ host, data = ini.tPn, FUN = function(x){paste0(x, collapse = ";")})[,2]
sim.df$transposon = apply(sim.df, 1, FUN = function(x){return(tPn.x(x[3]))})

# 3 duplicating between chromosomes
sim.df$transposon = apply(sim.df, 1, FUN = function(x){paste0(x[3], ";", tPn.flat(x[3], ori = 0))})
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
  save(hOst, rec.host, rec.transposon, rec.offspring, rec.gRecom, rec.tpnCount, sim.df, file = paste0("../data/tPn--", argv[3], "_", argv[5], ".rda"), compress = "xz")
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
  if(sum(s0.fit$fitness)==0){ eLim = -2; break }

  # 4 Chromosome segregation
  s0.fit = rbind(s0.fit, s0.fit)
  s0.fit$chromosome = rep(c(0,1), each = hPop)
  s0.tpnSeg = aggregate(tpn ~ host + chromosome, data = s0.tpn, FUN = function(x){return(paste(x, collapse = ";"))})
  s0.tpnSeg = merge(s0.tpnSeg, s0.fit, all = T)
  s0.tpnSeg = s0.tpnSeg[order(s0.tpnSeg$chromosome, s0.tpnSeg$host),]
  s0.tpnSeg$tpn[is.na(s0.tpnSeg$tpn)] = ""
  s0.fit$chromosome = s0.tpnSeg$fitness = NULL

  s0.tpnSeg$tpn = tPn.flat(s0.tpnSeg$tpn, ori = 1)
  s0.fit$host = s0.tpnSeg$host = seq_len(nrow(s0.fit))

  # 5 Repopulate
  fTree = matrix(rep(sample(s0.fit$host, hPop * ifelse(host.0$cell=="haploid", 1, 2), replace = T, prob = s0.fit$fitness/sum(s0.fit$fitness)), ifelse(host.0$cell=="haploid", 2, 1)), ncol = 2)

  # 6 Clean record dataframe
  s0 = sim.df
  for(i in seq_len(ncol(sim.df))){ sim.df[,i] = "" }

  # 7 Inherit appropriate parameters
  sim.df$familyTree = paste(fTree[,1], fTree[,2], sep = ";")
  sim.df$hARecomRate = h1.mod(rowMeans(matrix(s0$hARecomRate[match(fTree - (fTree > hPop) * hPop, as.numeric(row.names(s0)))], ncol = 2)), host.0$hAReomH1)
  tpns = matrix(s0.tpnSeg$tpn[match(fTree, s0.tpnSeg$host)], ncol = 2)
  tpns[,2] = tPn.flat(tpns[,2], ori = 0)
  sim.df$transposon = apply(tpns, 1, FUN = function(x){return(paste(x[nzchar(x)], collapse = ";"))})

  ##### Transposon activating stage #####
  cat(" Transposons;")
  # 1 Set transposon players
  s0.tpn = strsplit(gsub("^;$","", sim.df$transposon), ";")
  s0.tpnLen = lengths(s0.tpn)
  if(sum(s0.tpnLen)==0){ eLim = -1; break }
  s0.tpn = data.frame(tpn = unlist(s0.tpn), host = rep(seq_len(length(s0.tpn)), s0.tpnLen), fitness = 1)
  s0.tpn = sDF.replace = cbind(s0.tpn[s0.tpn$tpn!="",], tPn.io(s0.tpn$tpn))
  s0.tpn$numTpn = rep(s0.tpnLen, s0.tpnLen)
  s0.final = s0.tpn[!(s0.tpn$valid),] # spliced transposons are grounded
  s0.final$numTpn = NULL
  s0.tpn = s0.tpn[s0.tpn$valid,]

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

  # 3 Jump? Copy? Replicating number? (jump first, then copy)
  s0.tpn[,c("mv", "cp", "cpRep")] = matrix(c(runif(nrow(s0.tpn) * 2), rpois(nrow(s0.tpn), 1)), ncol = 3)
  while(any(s0.tpn$cpRep < 1)){s0.tpn$cpRep[which(s0.tpn$cpRep < 1)] = rpois(sum(s0.tpn$cpRep < 1), 1)}
  s0.tpn$mv = ifelse(s0.tpn$mv > s0.tpn$jumpRateFinal, NA, s0.tpn$mv)
  s0.tpn$cp = ifelse(s0.tpn$cp > s0.tpn$copyRateFinal, NA, s0.tpn$cp)

  # 4 Get transposons with activity
  s0.final = rbind(s0.final, s0.tpn[which(is.na(s0.tpn$mv)),colnames(s0.final)])
  tpns = c(which(!is.na(s0.tpn$mv)), rep(which(!is.na(s0.tpn$cp)), s0.tpn$cpRep[which(!is.na(s0.tpn$cp))]))
  if(length(tpns)>0){
    tpns = s0.tpn[tpns,]
    tpns = tpns[order(tpns$host),]
    tpns$replacement = NA

  # 5 Get replacement transposon for jumps
    if(any(!is.na(tpns$mv))){
      tpns$replacement = sim.df$hARecomRate[match(tpns$host, row.names(sim.df))]
      tpns$replacement = ifelse(is.na(tpns$mv), NA, tpns$replacement)
      tpns$replacement[which(!is.na(tpns$replacement))] = runif(sum(!is.na(tpns$replacement))) <= tpns$replacement[which(!is.na(tpns$replacement))]
      sim.df.fail = aggregate(replacement ~ host, data = tpns, FUN = prod)
      sim.tpn.fail = tpns[!is.na(tpns$replacement) & tpns$replacement==0,c("host", "gene")]
      if(nrow(sim.tpn.fail)>0){
        sim.tpn.fail = aggregate(gene ~ host, data = sim.tpn.fail, FUN = function(x){paste(x, collapse = ";")})
        sim.df$tpnFail = sim.tpn.fail$gene[match(as.numeric(row.names(sim.df)), sim.tpn.fail$host)]
        sim.df$tpnFail[is.na(sim.df$tpnFail)] = ""
      }
      ## Use homolog's transposons
      ok = unique(tpns[which(tpns$replacement==1), c("gene", "location", "host")])
      ok$repG = paste0(substr(ok$gene,1,1), (as.numeric(substr(ok$gene,2,2)) + 1) %% 2, substr(ok$gene,3,nchar(ok$gene)))
      hom = sDF.replace[paste(sDF.replace$host, sDF.replace$gene, sDF.replace$location) %in% paste(ok$host, ok$repG, ok$location),]
      if(nrow(hom)>0){
        cp = as.numeric(substr(hom$gene,2,2))
        hom$tpn = mapply(function(s, k){tPn.flat(s, ori = k)}, hom$tpn, cp, USE.NAMES = FALSE)
        hom$gene = paste0(substr(hom$gene,1,1), (cp + 1) %% 2, substr(hom$gene,3,nchar(hom$gene)))
        s0.final = rbind(s0.final, hom)
      }
    }

  # 6 Transposon relocation
    tpns$loc = runif(nrow(tpns))
    tpns$tpn = strsplit(tPn.io(tpns[,which(colnames(tpns) %in% colnames(tPn.0))]), ";")[[1]] # incorporating evolved jump/copy rates
    tpns = tPn.newloc(tpns, gen = gEn)

  ## 6.1 Copy direction allowed?
    if(any(tPn.0$copyDir!="both") & any(!is.na(tpns$cp))){
      tpns.cp = tpns[!is.na(tpns$cp),]
      tpns = tpns[is.na(tpns$cp),]
      tpns.cp$dIr = apply(tpns.cp, 1, FUN = function(x){return(tPn.comLoc(c(x[1], x[23])))}) == tpns.cp$copyDir
      tpns.cp$dIr[tpns.cp$copyDir=="both"] = TRUE
      tpns = rbind(tpns, tpns.cp[tpns.cp$dIr, which(colnames(tpns.cp) %in% colnames(tpns))]); rm(tpns.cp)
    }

  ## 6.2 New transposons inclusion
    tpns$tpn = tpns$newLoc
    tpns = cbind(tpns[,1:3], tPn.io(tpns$tpn))

  # 7 Harmonize all transposons from activity stage
    s0.final = rbind(s0.final, tpns[,which(colnames(tpns) %in% colnames(s0.final))])
  }
  if(nrow(s0.final)<1){ eLim = -1; break }

  ##### Gene recombination stage (assume no transposons excise / add complications) #####
  cat(" Recombination\n")
  # 1 Check if there is recombination happening
  reCom = TRUE
  reCom.df = as.data.frame(matrix(0, nrow = hPop, ncol = length(rG)))
  colnames(reCom.df) = rG
  for(i in seq_len(length(rG))){reCom.df[,i] = as.numeric(grepl(paste0("g0", colnames(reCom.df)[i]), sim.df$transposon, fixed = TRUE)) + as.numeric(grepl(paste0("g1", colnames(reCom.df)[i]), sim.df$transposon, fixed = TRUE))};rm(i)

  # 2 Number of genes going to be recombined
  gRecRate = host.0$recom * nrow(geneIn$gene) / 5000 # [202610231012] So 1/1000 genes undergo recombination at the highest rates; this is the product of the proportion of the population undergoing recombination, and the probability of a gene being affected by a recombination in such an isolate; this probability is the product of the number of recombinations, and the length of the recombinations (e.g. ~1% of pneumococci undergo transformation in culture; ~10% of a particular locus undergo recombination; if we estimate a mean length of 5kb, this implies ~40 recombinations per individual); so if you work out how many genes are affected by each recombination on average, you can scale this up such that an individual locus has a maximal recombination rate of ~0.001
  numGenes = lengths(strsplit(sim.df$transposon, ";"))
  nG.no = which(numGenes==0)
  numGenes = data.frame(tpn = tPn.io(paste(sim.df$transposon, collapse = ";"))$size, host = rep(seq_len(hPop), numGenes))
  numGenes = aggregate(tpn ~ host, data = numGenes, FUN = sum)
  if(length(nG.no)>0){
    nG.no = as.data.frame(matrix(c(nG.no, rep(0, length(nG.no))), ncol = 2))
    colnames(nG.no) = colnames(numGenes)
    numGenes = rbind(numGenes, nG.no)
    numGenes = numGenes[order(numGenes$host),]
  }
  numGenes = gRecRate * (mean(geneIn$gene$length) * (1 + numGenes$tpn / geneIn$genome[2] * (host.0$recomH1!="switch")) * as.numeric(!(apply(reCom.df, 1, max) == 2)))
  for(i in seq_len(hPop)){numGenes[i] = rpois(1, lambda = numGenes[i])}
  numGenes = pmin(numGenes, nrow(geneIn$gene))
  if(all(apply(reCom.df, 1, max) == 2) | gRecRate==0 | all(numGenes==0)){reCom = FALSE}

  if(reCom){
  # 3 Hotspots selection
  rGenes = c(); for(i in seq_len(hPop)){rGenes = c(rGenes, sample(geneIn$gene$locus_tag, numGenes[i], replace = FALSE))};rm(i)
  rGenes = data.frame(gene = rGenes, host = rep(seq_len(hPop), numGenes))
  hRecom = which(numGenes>0)

  # 4 Nearby gene regions
  rGSel = NULL; for(i in hRecom){ rGSel = rbind(rGSel, cbind(gRc.prob(rGenes$gene[rGenes$host==i]), i)) };rm(i)
  colnames(rGSel)[3] = "host"

  # 5 Probability of each gene in each host being selected
  rGenes = matrix(rep(geneIn$gene$length, hPop), ncol = hPop)
  tpns = aggregate(size ~ host + gene, data = s0.final, FUN = sum)
  tpns = tpns[substr(tpns$gene,1,1)=="g",]
  tpns$gene = match(substr(tpns$gene,2,nchar(tpns$gene)), geneIn$gene$locus_tag)
  for(i in seq_len(nrow(tpns))){rGenes[tpns$gene[i], tpns$host[i]] = rGenes[tpns$gene[i], tpns$host[i]] + tpns$size[i]};rm(i)
  rGenes = 1 / rGenes

  ## Remove broken genes in tpnFail
  if(any(sim.df$tpnFail!="")){
    s0.tpn = strsplit(sim.df$tpnFail, ";")
    for(i in which(lengths(s0.tpn)>0)){
      s0.tpn[[i]] = s0.tpn[[i]][substr(s0.tpn[[i]],1,1)=="g"]
      if(length(s0.tpn[[i]])>0){ rGenes[geneIn$gene$locus_tag %in% substr(s0.tpn[[i]],2,nchar(s0.tpn[[i]])),i] = 0 }
    };rm(i)
  }

  rGenes = rGenes / rowSums(rGenes)

  # 6 Gene recombination
  sDF.replace.success = s0.final[which(s0.final$host==0),] # initiate empty datafrmae inheriting structure of s0.final
  s0.donor = s0.final
  if(host.0$cell=="diploid"){
    rGSel$src = rGSel$host; rGSel$chromosome = (as.numeric(substr(rGSel$locus_tag, 1, 1)) + 1) %% 2
    rGSel$donor = paste0(rGSel$chromosome, substring(rGSel$locus_tag, 2))
  }else{
    rG.hap = unique(rGSel[,c("contiguous", "host")])
    rG.hap$src = NA
    rGSel$donor = rGSel$locus_tag
    for(i in seq_len(nrow(rG.hap))){
      rG.t = rGenes[which(geneIn$gene$locus_tag %in% rGSel$locus_tag[which(rGSel$contiguous==rG.hap$contiguous[i] & rGSel$host==rG.hap$host[i])]),-rG.hap$host[i], drop = FALSE]
      rG.hap$src[i] = sample(seq_len(hPop)[-rG.hap$host[i]], 1, prob = colSums(rG.t)/sum(rG.t)) # recombination gene source
    };rm(i,rG.t)
    rG.hap$chromosome = sample(c(0,1), nrow(rG.hap), replace = TRUE)
    rGSel$src = rG.hap$src[match(paste(rGSel$contiguous,rGSel$host), paste(rG.hap$contiguous,rG.hap$host))]
    rGSel$chromosome = rG.hap$chromosome[match(paste(rGSel$contiguous,rGSel$host), paste(rG.hap$contiguous,rG.hap$host))]
  }

  for(i in seq_len(nrow(rGSel))){
    g.tpn = paste0("g", rGSel$locus_tag[i])
    d.tpn = paste0("g", rGSel$donor[i])
    sim.df$tpnFail[rGSel$host[i]] = sub("^;", "", sub(";$", "", sub(g.tpn, "", sim.df$tpnFail[rGSel$host[i]]))) # tpnFail rescue

    # Transposon HGT
    sDF.replace = s0.donor[which(s0.donor$gene == d.tpn & s0.donor$host == rGSel$src[i]),]
    s0.final = s0.final[!(s0.final$gene == g.tpn & s0.final$host == rGSel$host[i]),]
    if(nrow(sDF.replace)>0){
      if(d.tpn != g.tpn){ sDF.replace$tpn = tPn.flat(sDF.replace$tpn, ori = rGSel$chromosome[i]); sDF.replace$gene = g.tpn }
      sDF.replace$host = rGSel$host[i]
      sDF.replace.success = rbind(sDF.replace.success, sDF.replace)
    }
  };rm(i, g.tpn, d.tpn)

  # 7 Intergenic recombination (if prior and post region gene recombination from the same source)
  rGSel = rGSel[order(rGSel$host, rGSel$chromosome, rGSel$locus_tag),]
  rGSel$gLoc = match(rGSel$locus_tag, geneIn$gene$locus_tag)
  rGSel$iReplace = FALSE
  for(i in seq_len(nrow(rGSel)-1)){
    if(rGSel$host[i]==rGSel$host[i+1] && rGSel$src[i]==rGSel$src[i+1] && rGSel$chromosome[i]==rGSel$chromosome[i+1] && rGSel$gLoc[i+1]==(rGSel$gLoc[i]+1)){rGSel$iReplace[i+1] = TRUE}
  };rm(i)
  if(any(rGSel$iReplace)){
    rGenes = rGSel[rGSel$iReplace,]
    rGenes$locus_tag = paste0("i", rGenes$locus_tag); rGenes$donor = paste0("i", rGenes$donor)
    for(i in seq_len(nrow(rGenes))){
      tmp = s0.donor[which(s0.donor$gene==rGenes$donor[i] & s0.donor$host==rGenes$src[i]),]
      s0.final = s0.final[!(s0.final$gene==rGenes$locus_tag[i] & s0.final$host==rGenes$host[i]),]
      if(nrow(tmp)>0){
        if(rGenes$donor[i] != rGenes$locus_tag[i]){
          tmp$tpn = tPn.flat(tmp$tpn, ori = rGenes$chromosome[i])
          tmp$gene = rGenes$locus_tag[i]
        }
        tmp$host = rGenes$host[i]; sDF.replace.success = rbind(sDF.replace.success, tmp)
      }
    };rm(i)
  }
  s0.final = rbind(s0.final, sDF.replace.success)
  if(nrow(s0.final)<1){ eLim = -1; break }

  # 8 Record gene recombinations (rec.host)
  rGSel$paste = paste(rGSel$gLoc, rGSel$src, rGSel$chromosome, sep = "!")
  rGSel = aggregate(paste ~ host, data = rGSel, FUN = function(x){paste(x, collapse = ";")})
  sim.df$host[hRecom] = rGSel$paste
  }

  ##### Harmonize all transposon records from all stages #####
  s0.final = aggregate(tpn ~ host, data = s0.final, FUN = function(x){paste0(x, collapse = ";")})
  sim.df$transposon = s0.final[match(as.numeric(row.names(sim.df)), s0.final$host),2]
  if(sum(is.na(sim.df$transposon))>0){sim.df$transposon[is.na(sim.df$transposon)] = ""}
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
