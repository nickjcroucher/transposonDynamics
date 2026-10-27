#!/bin/env Rscript
# author: ph-u
# script: ai_convergence.r
# desc: paths
# in: source("ai_convergence.r")
# out: NA
# arg: 0
# date: 20261005

source("src.r")

table(f.0[,-1])
f.ref[f.ref$Accession %in% f.0$Accession[f.0[,2]=="Transposase" & f.0[,3]=="Unknown"],]
f.ref[f.ref$Accession %in% f.0$Accession[f.0[,2]=="Unknown" & f.0[,3]=="Cellular functions"],]

f.1 = f.0[which(f.0[,2] %in% c("Reverse transcriptase", "Transposase") | f.0[,3] %in% c("Reverse transcriptase", "Transposase")),]
write.table(f.1, "../data/CL0219_category.csv", sep = ",", col.names = F, quote = F, row.names = F)

##### Env: CL0219 #####
d.0 = as.data.frame(fromJSON("../raw/CL0219_raw.jsonl"))
write.table(d.0$results.metadata$accession, "../raw/CL0219_tax_report_raw.txt", quote = F, row.names = F, col.names = F)
#d.0 = cbind(d.0$results.metadata, d.0$results.sets) # infinitely loading
colnames(d.0$results.metadata)
d.1 = as.data.frame(t(matrix(unlist(d.0$results.sets), nrow = 3)))
system("sed -e 's/|\t/@/g' -e 's/\t//g' -e 's/ @/@/g' ../raw/CL0219_tax_report.txt > ../raw/CL0219_tax_report-1.txt") # excel delimiter text to columns
d.2 = read.table("../raw/CL0219_tax_report-1.txt", header = T, sep = "\t", quote = "", comment.char = "")
write.table(unique(c(d.2$taxid, unlist(strsplit(d.2$lineage, " ")))), "../raw/CL0219_tax_report_rawFull.txt", quote = F, row.names = F, col.names = F)
system("sed -e 's/|\t/@/g' -e 's/\t//g' -e 's/ @/@/g' ../raw/CL0219_tax_report-Full.txt > ../raw/CL0219_tax_report-Full-1.txt") # excel delimiter text to columns
d.2 = read.table("../raw/CL0219_tax_report-Full-1.txt", header = T, sep = "@", quote = "", comment.char = "", fill = TRUE)
d.2.lineage = strsplit(d.2$lineage, " ")
d.2.df = as.data.frame(matrix(nrow = nrow(d.2), ncol = max(lengths(d.2.lineage))))
cat(as.character(t.ori <- Sys.time()), ": CL0219 domain distribution lineage collection\n"); for(i in seq_len(nrow(d.2))){
  eTa(t.ori, i, nrow(d.2))
  d.2.df[i,seq_len(length(d.2.lineage[[i]]))] = rev(d.2.lineage[[i]])
};rm(i);cat("\n") # 7 mins

##### Load CL0219 taxonomic distributions #####
f.json = f.all[grep("/CL0219/", f.all)]
f.pfam = strsplit(basename(f.json), "-")
cat(as.character(t.ori <- Sys.time()), ": Loading CL0219 domain distributions\n"); for(i in seq_len(length(f.json))){
  eTa(t.ori, i, length(f.json))
  i0 = flat(read_json(f.json[i])$taxa)
  i0 = i0[which(i0$rank=="domain"),]
  i0$file = f.pfam[[i]][1]
  if(i>1){d.0 = rbind(d.0, i0)}else{d.0 = i0}
};rm(i, i0);cat("\n")
# i1 = 0; i0 = c(); for(i in seq_len(length(f.json))){
#   cat(date(), ":", f.pfam[[i]][1], "\n")
#   i0 = try(fromJSON(f.json[i], flatten = TRUE, simplifyDataFrame = TRUE))
#   if(class(i0)!="try-error"){ #!!!
#     i2 = as.data.frame(matrix("", nrow = length(i0$taxa), ncol = unique(lengths(i0$results$entries))))
#     colnames(i2) = names(i0$results$entries[[1]])
#     for(i3 in seq_len(length(i0$results$entries))){ i2[i3,] = i0$results$entries[[i3]] }
#     i2 = cbind(i0$results$metadata, i2)
#     i2$file = f.pfam[[i]][1]
#     i1 = i1 + 1
#     if(i1>1){d.0 = rbind(d.0, i2)}else{d.0 = i2}
#   }else{print(f.pfam[[i]][1])}
# };rm(i, i0, i1, i2, i3)

##### Categorize protein domains in taxonomy #####
d.0$category = "Prokaryotes"
d.0$category[d.0$name %in% c("Eukaryota")] = "Eukaryotes"
d.1 = d.0[which(d.0$name %in% c("Archaea", "Bacteria", "Eukaryota")), c(5:9)]
# d.0$category = d.2.df[match(d.0$accession, d.2$taxid),1]
# d.0$category[is.na(d.0$category)] = -1
# d.0$category[d.0$category == 131567] = d.2.df[match(d.0$accession[d.0$category == 131567], d.2$taxid),2]
# d.0$category = d.2$taxname[match(d.0$category, d.2$taxid)]
# d.0$category[d.0$category %in% c("Archaea", "Bacteria")] = "Prokaryotes"
d.1$type = f.1$Classification.x[match(d.1$file, f.1$Accession)]

##### Correlation #####
# x0 = table(d.0[d.0$category %in% c("Prokaryotes", "Eukaryota"),c("category", "type")])
d.2 = aggregate(proteins ~ category + type, data = d.1, FUN = sum)
x0 = matrix(d.2[,3], nrow = length(unique(d.2$category))); row.names(x0) = unique(d.2$category); colnames(x0) = unique(d.2$type)
x0;chisq.test(x0)
d.2 = aggregate(species ~ category + type, data = d.1, FUN = sum)
x0 = matrix(d.2[,3], nrow = length(unique(d.2$category))); row.names(x0) = unique(d.2$category); colnames(x0) = unique(d.2$type)
x0;chisq.test(x0)
