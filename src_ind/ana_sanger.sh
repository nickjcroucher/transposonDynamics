#!/bin/env bash
# author: ph-u
# script: ana_sanger.sh
# desc: HPC script for ana_ind_child.sh
# in: bash ana_sanger.sh
# out: NA
# arg: 0
# date: 20260702

mkdir -p ../data

sed -e "s/mAx/$(( `wc -l < ../upload/scenario.csv` - 1 ))/" ana_ind_child.sh > sC.sh

bsub < sC.sh

exit
