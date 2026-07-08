#!/bin/env bash
# author: ph-u
# script: sangerHPC.sh
# desc: HPC script for transposon model
# in: bsub sangerHPC.sh
# out: NA
# arg: 0
# date: 20260702

PATH="/software/isg/private/wrappers/apptainer/1.4.0:$PATH"

mkdir -p ../data && apptainer pull docker://ghcr.io/nickjcroucher/transposondynamics:latest

bsub < sangerHPC_child.sh

exit
