#!/bin/env bash
# author: ph-u
# script: setup.sh
# desc: setup openstack instance
# in: bash setup.sh
# out: directory setup
# arg: 0
# date: 20261027

mkdir -p pj02/{data,src,work}

cd pj02/data
rsync -avH ph17@farm22.internal.sanger.ac.uk:/data/pam/team377/ph17/scratch/pj02/openstack/ .
mv *.sh ../src

cd ../src
bash sangerHPC.sh

exit
