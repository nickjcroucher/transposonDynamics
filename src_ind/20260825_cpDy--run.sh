#!/bin/env bash
# author: ph-u
# script: 20260825_cpDy--run.sh
# desc: set script ro run Rmd rendering
# in: bash 20260825_cpDy--run.sh
# out: NA
# arg: 0
# date: 20260826

Rscript 20260825_cpDy--list.r
sed -e "s/mAx/`wc -l < ../upload/20260825_cpDy--list.csv`/" 20260825_cpDy.sh > sC.sh

bsub < sC.sh

exit
