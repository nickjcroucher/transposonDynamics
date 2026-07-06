#!/bin/env bash
# author: ph-u
# script: simMaster.sh
# desc: run simulation
# in: bash simMaster.sh [../raw/input.csv] [../raw/hpcTag.csv] [../raw/seed.csv] [../raw/scenario.csv] [hpcTag line number]

# for i in `seq 1 36`;do bash simMaster.sh ../raw/input.csv ../raw/hpcTag.csv ../raw/seed.csv ../raw/scenario.csv ${i} > ${i}.e.txt & done
# for i in `ls *.e.*`;do [[ `grep -e "warnings" ${i}` ]] && echo -e "${i}";done

[[ -z $5 ]] && echo -e "Need 5 inputs\n`head -n 5 $0 | tail -n 1`" && exit

inpCSV=$1
hpcCSV=$2
sedCSV=$3
sceCSV=$4
hPc=$(( $5 + 1 ))

hpcTag=`head -n ${hPc} ${hpcCSV} | tail -n 1`
sEed=`echo -e "${hpcTag}" | cut -f 1 -d ","`
sCen=`echo -e "${hpcTag}" | cut -f 2 -d ","`

Rscript simulate.r ${inpCSV} ${sedCSV} ${sEed} ${sceCSV} ${sCen}

exit
