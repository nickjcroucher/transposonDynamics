#!/bin/env Rscript
# author: ph-u
# script: pfam_analysis.r
# desc: pfam CL0219 domain analysis
# in: Rscript pfam_analysis.r
# out: <stdout>
# arg: 0
# date: 20261123

##### data import #####
cat(date(),": Data import\n")
f0 = list.files("../CL0219-pfam/scripts/tables/", full.names = TRUE)
d.clan_families = read.table(f0[grep("clan_families", f0)], sep = "\t") # 8.6 KB
d.proteome_citations = read.table(f0[grep("proteome_citations", f0)], sep = "\t", comment.char = "", quote = "") # 36 MB
d.proteomes = read.table(f0[grep("proteomes", f0)], sep = "\t", comment.char = "", quote = "") # 153 MB
d.taxonomy_lineage = read.table(f0[grep("taxonomy_lineage", f0)], sep = "\t", comment.char = "", quote = "") # 36 MB
d.matches = read.delim(gzfile(f0[grep("matches", f0)])) # 14 MB
d.protein_proteome = read.delim(gzfile(f0[grep("protein_proteome", f0)])) # 9.8 MB

##### Analysis #####
cat(date(),": Done\n")
