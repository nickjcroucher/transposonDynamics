#!/bin/env bash
# author: ph-u
# script: ref_genome_prep.sh
# desc: download genome from ENA downloader and format genome length data
# in: bash ref_genome_prep.sh [../custom/loc/input.csv]
# out: data/ref_genome_prep.csv
# arg: 1 (optional)
# date: 20260622

[ -z $1 ] && inFile="../raw/input.csv" || inFile=$1

aCc=`head -n 2 ${inFile} | tail -n 1 | cut -f 2 -d "," | cut -f 1 -d "."`

##### Genome Download #####

mkdir -p ../data

#java -jar ena-file-downloader.jar --accessions=${aCc} --format=READS_FASTQ --location=../data/ --protocol=FTP --asperaLocation=null # https://github.com/natesheehan/rena/tree/main/inst/extdata/ena-file-downloader

[ -d ../data/ref ] && rm -r ../data/ref/

datasets download genome accession ${aCc} --dehydrated --filename ref.zip --include gff3 && unzip ref.zip -d ../data/ref && datasets rehydrate --directory ../data/ref/

##### Gene annotation file extraction #####

rEf=`find ../data/ref/ | grep -e ".gff$"`
echo -e "${rEf}"
mv ${rEf} ../data/${aCc}.gff

rm ref.zip
rm -r ../data/ref

exit
