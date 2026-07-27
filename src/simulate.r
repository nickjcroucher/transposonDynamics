#!/bin/env Rscript
# author: ph-u
# script: simulate.r
# desc: Discrete-time model for transposon ecology within one genome
# in: Rscript simulate.r [../custom/loc/input.csv] [../custom/loc/seed.csv] [seed number] [../custom/loc/scenario.csv] [scenario number]
# out: NA
# arg: 1 (optional)
# date: 20260623

argv=(commandArgs(T))
if(length(argv) < 5){argv=c("../raw/input.csv", "../raw/seed.csv", "1", "../raw/scenario.csv", "7")}

##### env set #####
cat(date(),": set environment",argv[3],"-",argv[5],"\n")
source("func.r")
set.seed(read.csv(argv[2], header = F)[,1][as.numeric(argv[3])])
inFile = inParams(argv[1])
sCene = read.csv(argv[4], header = T)[as.numeric(argv[5]),]
gPrm = c(toxicProb = as.numeric(inFile$params$Value[inFile$params$Type=="percentage transposon perturbation genotoxic"])/100,
           boostProb = as.numeric(inFile$params$Value[inFile$params$Type=="percentage chance transposon perturbation boost"])/100,
           boostCoef = as.numeric(inFile$params$Value[inFile$params$Type=="percentage amplitude transposon perturbation boost"])/100)

##### Initiate populations #####
cat(date(),": initiate population",argv[3],"-",argv[5],"\n")
hOst = ini.host(inFile$params$Value[inFile$params$Type=="host genome variation"], inFile$gene, inFile$params$Value[inFile$params$Type=="host genetic variation"])
tPn = ini.transposon(inFile$params$Value[inFile$params$Type=="transposon size in bp"], scenario = sCene)

##### Initiate record dataframes #####
gEn.max = as.numeric(inFile$params$Value[inFile$params$Type=="host organism constant generation number"])

perturbGen = ceiling(rev(seq(1, gEn.max, gEn.max/(sCene$genotoxic+1))))-1
perturbGen[perturbGen==0] = gEn.max + 1 # a numeric placeholder that can never achieve

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
    g.tmp = tPn.act(tPn = tPn.tag[i0], gen = 0, gene.df = g.tmp, scenario = sCene, pAram = gPrm, gToxic = F)
    tPn.loc[i0] = tPn.get(g.tmp)
    colnames(g.tmp)[1] = tPn.get(g.tmp, F)
  };rm(i0, g.tmp)

  ## Validate transposons
  if(length(tPn.loc) > 1){for(i0 in 1:(length(tPn.loc)-1)){for(i1 in (i0+1):length(tPn.loc)){
    tPn.loc[i0] = tPn.x(tPn1 = tPn.loc[i0], tPn2 = tPn.loc[i1])
  }};rm(i0,i1)}
  sim.df$transposon[i] = paste0(tPn.loc, collapse = ";")
};rm(i, tPn.loc)

##### Wright-Fisher / Neutral model run #####
tPn.iDx = "transposon population size per genome sd"
gEn = 0; repeat{
  if((gEn %% 10) == 0){cat(date(),": run simulation",argv[3],"-",argv[5], "; generation =",gEn, "/",gEn.max, "(",round(gEn/gEn.max*100,2),"% )\n")}
  ## Population dynamics snapshot
  rec.host[gEn + 1,] = sim.df$host
  rec.transposon[gEn + 1,] = sim.df$transposon
  rec.offspring[gEn + 1,] = sim.df$familyTree

  gEn = gEn + 1; if(gEn > gEn.max){ break } # simulation done
#  cat(date(),": generation",gEn,"\n")

  ## Population reproduction stage
  sim.df = host.reproduce(res.pool = sim.df, gene.df = inFile$gene, transposon.size = inFile$params$Value[inFile$params$Type=="transposon size in bp"], fitness.advantage = as.numeric(inFile$params$Value[inFile$params$Type=="percentage of fitness benefit with transposon"]))

  ## Transposon jumping stage
  ### 1. Set jumping indicators
  tPn.sums = lengths(tPn.list <- strsplit(sim.df$transposon, ";"))
  tPn.csum = cumsum(tPn.sums)
  for(i in 1:nrow(sim.df)){

  ### 2. Reconstruct transposon-inserted gene table
    g.tmp = reGeneDF(
      tPn = sim.df$transposon[i],
      tPn.size = as.numeric(inFile$params$Value[inFile$params$Type=="transposon size in bp"]),
      gene.df = inFile$gene
      )

    if(tPn.sums[i] > 0){
  ### 3. Map indicators with transposon locations
      tPn.tag = tPn.io(sim.df$transposon[i])
      for(i1 in 1:tPn.sums[i]){
        if(class(tPn.tag)=="character"){
          g.tmp = tPn.act(tPn = tPn.io(tPn.tag), gen = gEn, gene.df = g.tmp, scenario = sCene, pAram = gPrm, gToxic = (gEn %in% perturbGen))
        }else{
          g.tmp = tPn.act(tPn = tPn.io(tPn.tag[i1,]), gen = gEn, gene.df = g.tmp, scenario = sCene, pAram = gPrm, gToxic = (gEn %in% perturbGen))
        }
        tPn.list[[i]][i1] = tPn.get(g.tmp)
        colnames(g.tmp)[1] = tPn.get(g.tmp, F)
      };rm(i1)

  ### 4. Validate each transposon
      tPn.list[[i]] = tPn.r(tPn.list[[i]])
      if(length(tPn.list[[i]])>1){  for(i1 in 1:(length(tPn.list[[i]])-1)){ for(i2 in (i1+1):length(tPn.list[[i]])){
          tPn.list[[i]][i1] = tPn.x(tPn1 = tPn.list[[i]][i1], tPn2 = tPn.list[[i]][i2])
      }};rm(i1,i2) }
      # tPn.list[[i]] = tPn.r(paste0(tPn.list[[i]], collapse = ";"))
      sim.df$transposon[i] = paste0(tPn.list[[i]], collapse = ";")
    }
  };rm(i)

  ## Gene recombination stage (assume no transposons excise / add complications)
  sim.df = g.Recom(
    res.pool = sim.df,
    gene.df = inFile$gene,
    recomRate = sCene$gene[1]
    )

}

cat(date(),": printing warnings",argv[3],"-",argv[5],"\n")
print(warnings())

##### Simulation record export #####
cat(date(),": result export",argv[3],"-",argv[5],"\n")
save(rec.host, rec.transposon, rec.offspring, file = paste0("../data/tPn--", argv[3], "_", argv[5], ".rda"), compress = "xz")
cat(date(),": simulation completed",argv[3],"-",argv[5],"\n")
