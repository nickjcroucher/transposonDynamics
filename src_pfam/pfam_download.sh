#!/bin/en bash
# author: ph-u
# script: pfam_download.sh
# desc: download pfam InterPro taxonomy distribution
# in: bash pfam_download.sh
# out: CL0219/[pfam]-tax.json
# arg: 0
# date: 20261019

i0="CL0219"
mkdir -p ../data/${i0}
while read -r l;do
  i=`echo -e "${l}" | cut -f 1 -d ","`
  curl "https://www.ebi.ac.uk/interpro/api/entry/pfam/${i}/?taxa" > ../data/${i0}/${i}-tax.json &
#  python3 script-interPro.py ${i} > ../data/${i0}/${i}-tax.json & # data too heavy
done < ../data/CL0219_category.csv

exit

