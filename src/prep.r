#argv = (commandArgs(T))
argv = "../raw/accession.txt"

##### packages #####
pKgs = c("devtools")
iNstall = pKgs[!(pKgs %in% row.names(installed.packages()))]
pKgs = pKgs[pKgs %in% row.names(installed.packages())]
for(i in iNstall){
  install.packages(i, dependencies = T)
  library(i, character.only = T)
};rm(i, iNstall)
for(i in pKgs){library(i, character.only = T)};rm(i, pKgs)

##### base sequence download #####
dNa = read.table(argv[1])
