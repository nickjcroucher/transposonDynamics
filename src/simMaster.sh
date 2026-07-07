#!/bin/env bash
# author: ph-u
# script: simMaster.sh
# desc: run simulation with replicates
# in: bash simMaster.sh [../raw/input.csv] [../raw/seed.csv] [../raw/scenario.csv] [scenario line number]

# for i in `seq 1 36`;do bash simMaster.sh ../raw/input.csv ../raw/seed.csv ../raw/scenario.csv ${i} > ${i}.e.txt & done
# for i in `ls *.e.*`;do [[ `grep -e "warnings" ${i}` ]] && echo -e "${i}";done

[[ -z $4 ]] && echo -e "Need 4 inputs\n`head -n 5 $0 | tail -n 1`" && exit

inpCSV=$1
sedCSV=$2
sceCSV=$3
hPc=$4

for i in `seq 1 10`;do
  Rscript simulate.r ${inpCSV} ${sedCSV} ${i} ${sceCSV} ${hPc}
done

exit
