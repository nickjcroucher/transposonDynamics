#!/bin/en bash
# author: ph-u
# script: pfam_download.sh
# desc: download pfam InterPro taxonomy distribution
# in: bash pfam_download.sh
# out: CL0219/[pfam]-tax.json
# arg: 0
# date: 20261019

i0="CL0219"
mkdir -p ${i0}
while read -r l;do
  i=`echo -e "${l}" | cut -f 1 -d ","`
  python3 script-interPro.py ${i} > ${i0}/${i}-tax.json &
done < CL0219_category.csv

exit

