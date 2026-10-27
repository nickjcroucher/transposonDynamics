#!/bin/env bash
# author: ph-u
# script: sangerHPC.sh
# desc: HPC script for transposon model
# in: bsub sangerHPC.sh
# out: NA
# arg: 0
# date: 20260702

PATH="/software/isg/private/wrappers/apptainer/1.4.0:$PATH"

mkdir -p ../{data,work}

[[ -f transposondynamics_latest.sif ]] && rm transposondynamics_latest.sif
apptainer pull docker://ghcr.io/nickjcroucher/transposondynamics:latest

if [[ `command -v bsub >/dev/null 2>&1` ]];then
  sed -e "s/mAx/$(( `wc -l < ../data/rerun.csv` ))/" sangerHPC_child.sh > sC.sh
  bsub < sC.sh
else
  sed -e "s/LSB_JOBINDEX/1/" sangerHPC_child.sh > sC.sh
  i0=`wc -l < ../data/rerun.csv`
  for i in `seq 1 ${i0}`;do
    bash sC.sh ${i} 1> ../work/pj02-${i}.txt &
  done
fi

exit
