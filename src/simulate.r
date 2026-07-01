#!/bin/env Rscript
# author: ph-u
# script: simulate.r
# desc: Discrete-time model for transposon ecology within one genome
# in: Rscript simulate.r [../custom/loc/input.csv]
# out: NA
# arg: 1 (optional)
# date: 20260623

argv=(commandArgs(T))
if(length(argv) < 5){argv=c("../raw/input.csv", "../raw/seed.csv", "1", "../raw/scenario.csv", "1")}

##### env set #####
cat(date(),": set environment\n")
source("func.r")
set.seed(read.csv(argv[2], header = F)[,1][as.numeric(argv[3])])
inFile = inParams(argv[1])
sCene = read.csv(argv[4], header = T)[as.numeric(argv[5]),]

##### Initiate populations #####
cat(date(),": initiate population\n")
hOst = ini.host(inFile$params$Value[inFile$params$Type=="host genome variation"], inFile$gene, inFile$params$Value[inFile$params$Type=="host genetic variation"])
tPn = ini.transposon(inFile$params$Value[inFile$params$Type=="transposon size in bp"])

##### Initiate record dataframes #####
gEn.max = as.numeric(inFile$params$Value[inFile$params$Type=="host organism constant generation number"])
rec.host = rec.transposon = rec.offspring = as.data.frame(matrix(NA, nrow = gEn.max + 1, ncol = as.numeric(inFile$params$Value[inFile$params$Type=="host organism constant population size"])))
sim.df = as.data.frame(matrix(NA, nrow = ncol(rec.host), ncol = 3))
colnames(sim.df) = c("host", "transposon", "familyTree")

##### Initial population #####
sim.df$familyTree = sample(1:length(hOst), nrow(sim.df), replace = T)
sim.df$host = hOst[sim.df$familyTree]
for(i in 1:nrow(sim.df)){
  tPn.loc = lOc = rNumVec(f = "uniform", L = inFile$transposon.titre[i], p1 = 0, p2 = 1)
  tPn.tag = sample(tPn$ini, length(lOc), replace = T)

  ## Map transposons
  for(i0 in 1:length(lOc)){
    if(i0==1){g.tmp = as.data.frame(inFile$gene)}
    g.tmp = tPn.jump(tPn.tag = tPn.tag[i0], tPn.prob = lOc[i0], tPn.size = tPn$size[match(tPn.tag[i0], tPn$ini)], gene.df = g.tmp, tPn.move = 0, jumProb = 1, generation = 0)
    tPn.loc[i0] = strsplit(colnames(g.tmp)[1], "[.]")[[1]][1]
    colnames(g.tmp)[1] = strsplit(colnames(g.tmp)[1], "[.]")[[1]][2]
  };rm(i0, g.tmp)

  ## Validate transposons
  if(length(tPn.loc) > 1){for(i0 in 1:(length(tPn.loc)-1)){for(i1 in 2:length(tPn.loc)){
    tPn.loc[i0] = tPn.x(
      tPn1 = tPn.loc[i0],
      tPn1.len = as.numeric(inFile$params$Value[inFile$params$Type=="transposon size in bp"]),
      tPn2 = tPn.loc[i1],
      tPn2.len = as.numeric(inFile$params$Value[inFile$params$Type=="transposon size in bp"])
    )
  }};rm(i0,i1)}
  sim.df$transposon[i] = paste0(tPn.loc, collapse = ";")
};rm(i, tPn.loc)

##### Wright-Fisher / Neutral model run #####
cat(date(),": run simulation\n")
gEn = 0; repeat{
  ## Population dynamics snapshot
  rec.host[gEn + 1,] = sim.df$host
  rec.transposon[gEn + 1,] = sim.df$transposon
  rec.offspring[gEn + 1,] = sim.df$familyTree

  gEn = gEn + 1; if(gEn > gEn.max){ break } # simulation done
  cat(date(),": generation",gEn,"\n")

  ## Population reproduction stage
  sim.df = host.reproduce(res.pool = sim.df, gene.df = inFile$gene, transposon.size = inFile$params$Value[inFile$params$Type=="transposon size in bp"])

  ## Transposon jumping stage
  ### 1. Set jumping indicators
  tPn.iDx = "transposon population size per genome sd"
  tPn.sums = lengths(tPn.list <- strsplit(sim.df$transposon, ";"))
  tPn.csum = cumsum(tPn.sums)
  tPn.pRob = data.frame(
    jump = rNumVec(
      f = inFile$params$Value[inFile$params$Type==sub(" sd", " distribution", tPn.iDx)],
      L = max(tPn.csum),
      p1 = inFile$params$Value[inFile$params$Type==sub(" sd", " mean", tPn.iDx)],
      p2 = inFile$params$Value[inFile$params$Type==tPn.iDx]
      ),
    land = rNumVec(f = "uniform", L = max(tPn.csum), p1 = 0, p2 = 1))
  for(i in 1:nrow(sim.df)){

  ### 2. Reconstruct transposon-inserted gene table
    g.tmp = reGeneDF(
      tPn = sim.df$transposon[i],
      tPn.size = as.numeric(inFile$params$Value[inFile$params$Type=="transposon size in bp"]),
      gene.df = inFile$gene
      )

    if(tPn.sums[i] > 0){
  ### 3. Map indicators with transposon locations
      i0 = tPn.pRob[(ifelse(i>1,tPn.csum[i-1],0)+1):tPn.csum[i],]
      for(i1 in 1:nrow(i0)){
        g.tmp = tPn.jump(
          tPn.tag = tPn.list[[i]][i1],
          tPn.prob = i0$land[i1],
          tPn.size = as.numeric(inFile$params$Value[inFile$params$Type=="transposon size in bp"]),
          gene.df = g.tmp,
          tPn.move = i0$jump[i1],
          jumProb = sCene$transposon[1],
          generation = gEn
        )
        tPn.list[[i]][i1] = strsplit(colnames(g.tmp)[1], "[.]")[[1]][1]
        colnames(g.tmp)[1] = strsplit(colnames(g.tmp)[1], "[.]")[[1]][2]
      }

  ### 4. Validate each transposon
      tPn.list[[i]] = tPn.r(tPn.list[[i]])
      if(length(tPn.list[[i]])>1){ for(i1 in 1:(length(tPn.list[[i]])-1)){ for(i2 in 2:length(tPn.list[[i]])){
        tPn.list[[i]][i1] = tPn.x(
          tPn1 = tPn.list[[i]][i1],
          tPn1.len = as.numeric(inFile$params$Value[inFile$params$Type=="transposon size in bp"]),
          tPn2 = tPn.list[[i]][i2],
          tPn2.len = as.numeric(inFile$params$Value[inFile$params$Type=="transposon size in bp"])
          )
      }}}
      sim.df$transposon[i] = paste0(tPn.list[[i]], collapse = ";")
    }
  };rm(i, i0, i1)

  ## Gene recombination stage (assume no transposons excise / add complications)
  sim.df = g.Recom(
    res.pool = sim.df,
    gene.df = inFile$gene,
    recomRate = sCene$gene[1]
    )

}

##### Simulation record export #####
cat(date(),": result export\n")
save(rec.host, rec.transposon, rec.offspring, file = paste0("../data/tPn--", argv[3], "_", argv[5], ".rda"), compress = "xz")
cat(date(),": simulation completed\n")
