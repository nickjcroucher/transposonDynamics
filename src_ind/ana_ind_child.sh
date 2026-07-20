#!/bin/env bash
# author: ph-u
# script: ana_ind_child.sh
# desc: HPC script for summarizing data from the transposon model
# in: bsub ana_ind_child.sh
# out: NA
# arg: 0
# date: 20260716

# BSUB -G team377f
# BSUB -o ../work/pj02-%J-%I.o
# BSUB -e ../work/pj02-%J-%I.e
# BSUB -q week
# BSUB -M 2000
# BSUB -R "select[mem>2000] rusage[mem=2000] span[hosts=1]"
# BSUB -J "ind[1-720]"

Rscript ana_individual.r ${LSB_JOBINDEX}

exit
## any BSUB -J with the correct format, no matter at any line, later line will replace earlier lines # BSUB -J "tdy[1-360]"
