#!/bin/env bash
# author: ph-u
# script: 20260825_cpDy.sh
# desc: Rmd rendering different versions
# in: bsub < 20260825_cpDy.sh
# out: upload/20260825_cpDy--*.html
# arg: 0
# date: 20260826

#BSUB -G team377f
#BSUB -o ../work/pj02-%J-%I.o
#BSUB -e ../work/pj02-%J-%I.e
#BSUB -q normal
#BSUB -M 16000
#BSUB -R "select[mem>16000] rusage[mem=16000] span[hosts=1]"
#BSUB -J "cpDy[1-mAx]%10"

##### Setup ImageMagick #####
mkdir -p ~/.imagemagick
cat > ~/.imagemagick/policy.xml <<'EOF'
<policymap>
  <policy domain="resource" name="list-length" value="100000"/>
  <policy domain="resource" name="memory"      value="4GiB"/>
  <policy domain="resource" name="map"         value="4GiB"/>
  <policy domain="resource" name="disk"        value="50GiB"/>
  <policy domain="resource" name="file"        value="4096"/>
  <policy domain="resource" name="thread"      value="1"/>
</policymap>
EOF

export MAGICK_CONFIGURE_PATH=$HOME/.imagemagick
export MAGICK_THREAD_LIMIT=1
export MAGICK_TMPDIR=/data/pam/team377/ph17/scratch/pj02/.imagemagickTMP/${LSB_JOBINDEX}
mkdir -p "$MAGICK_TMPDIR"
identify -list resource        # confirm your values took effect
trap 'rm -rf "$MAGICK_TMPDIR"' EXIT

##### Rmd render #####
i=`head -n ${LSB_JOBINDEX} ../upload/20260825_cpDy--list.csv | tail -n 1`
p1=`echo -e ${i} | cut -f 1 -d ","`
p2=`echo -e ${i} | cut -f 2 -d ","`
p3=`echo -e ${i} | cut -f 3 -d ","`

sed -e "s/pAra1/${p1}/" -e "s/pAra2/${p2}/" -e "s/pAra3/${p3}/" 20260825_cpDy.Rmd > 20260825_cpDy--${p1}_${p2}_${p3}.Rmd

date
Rscript -e "options(gganimate.nframes = 501, gganimate.dev_args = list(width = 500, height = 350, units = 'px', res = 96)); rmarkdown::render('20260825_cpDy--${p1}_${p2}_${p3}.Rmd', output_dir='../upload', intermediates_dir='../work')"
date

rm 20260825_cpDy--${p1}_${p2}_${p3}.Rmd

exit

