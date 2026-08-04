#!/bin/env bash
# author: ph-u
# script: sangerHPC.sh
# desc: HPC script for transposon model
# in: bsub sangerHPC.sh
# out: NA
# arg: 0
# date: 20260702

PATH="/software/isg/private/wrappers/apptainer/1.4.0:$PATH"

mkdir -p ../data

[[ -f transposondynamics_latest.sif ]] && rm transposondynamics_latest.sif && apptainer pull docker://ghcr.io/nickjcroucher/transposondynamics:latest

sed -e "s/mAx/$(( `wc -l < ../data/rerun.csv` ))/" sangerHPC_child.sh > sC.sh

bsub < sC.sh

exit
