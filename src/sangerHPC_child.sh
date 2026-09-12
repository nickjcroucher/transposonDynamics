#!/bin/env bash
# author: ph-u
# script: sangerHPC_child.sh
# desc: HPC script for transposon model
# in: bsub sangerHPC_child.sh
# out: NA
# arg: 0
# date: 20260702

# BSUB -G team377f
# BSUB -o ../work/pj02-%J-%I.o
# BSUB -e ../work/pj02-%J-%I.e
# BSUB -q week
# BSUB -M 10000
# BSUB -R "select[mem>10000] rusage[mem=10000] span[hosts=1]"
# BSUB -J "tdy[1-mAx]"

PATH="/software/isg/private/wrappers/apptainer/1.4.0:$PATH"

apptainer run --bind ${PWD}/../data:/data --pwd /src transposondynamics_latest.sif bash simMaster.sh ../raw/input.csv ../raw/seed.csv ../raw/scenario.csv ${LSB_JOBINDEX}

exit
## any BSUB -J with the correct format, no matter at any line, later line will replace earlier lines # BSUB -J "tdy[1-360]"
# interactive session: apptainer shell --nv --env PREPEND_PATH=/opt/conda/bin --env R_LIBS_USER=/scratch/$USER/Rlib --bind ${PWD}/../data:/data --pwd /src transposondynamics_latest.sif "$@"
