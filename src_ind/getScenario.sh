#!/bin/env bash
# author: ph-u
# script: getScenario.sh
# desc: get filenames of pdfs to open
# in: bash getScenario.sh [jp] [cp] [rc]
# out: stdout
# arg: 3
# date: 20261103

i0=`Rscript getScenario.r $1 $2 $3 | tail -n 1`
echo -e "${i0}"
for i in ${i0};do
# open /Users/ph17/0_lustre/pj02/data/ana--1-${i}.pdf
 open /Users/ph17/0_lustre/pj02/os_data/ana--1-${i}.pdf
done

exit

