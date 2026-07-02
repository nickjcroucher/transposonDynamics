#!/bin/env bash
# author: ph-u
# script: sangerHPC.sh
# desc: HPC script for transposon model
# in: bsub sangerHPC.sh
# out: NA
# arg: 0
# date: 20260702

# BSUB -G team377f
# BSUB -o ../work/pj02-%J-%I.o
# BSUB -e ../work/pj02-%J-%I.e
# BSUB -q week
# BSUB -M 4000
# BSUB -R "select[mem>4000] rusage[mem=4000] span[hosts=1]"
# BSUB -J "tdy[1-360]"

PATH="/software/isg/private/wrappers/apptainer/1.4.0:$PATH"

apptainer pull tdy.sif oras://ghcr.io/nickjcroucher/transposondynamics:latest

apptainer exec tdy.sif bash simMaster.sh ../raw/input.csv ../raw/hpcTag.csv ../raw/seed.csv ../raw/scenario.csv ${LSB_JOBINDEX}

exit
