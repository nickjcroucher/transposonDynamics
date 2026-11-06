#!/bin/env Rscript
# author: ph-u
# script: getScenario.r
# desc: Extract scenario code from list
# in: Rscript getScenario.r [jp] [cp] [rc]
# out: stdout
# arg: 3
# date: 20261103

argv=(commandArgs(T))
source("../src/setScenario.r")
argv = as.numeric(argv)
cat(unique(sceSet(jp=argv[1], cp=argv[2], rc=argv[3], cT=unique(a0$cell), hAR=.5)[,2]), "\n")
