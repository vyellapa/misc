#!/bin/bash

DATE=$(date -d '1 hour ago' '+%Y-%m-%d')
CORES=4

echo "retrieving the bams..."
allbam=$(lk get_data -fi projects 222 -fi hasbam True)
#allbam="/ifs/res/leukgen/local/opt/leukdc/data/workflows/44/16/14416/data/bam/E-H-110915-T1-1-D1-1.bam"

mkdir -p cnvkit_${DATE}
mkdir -p cnvkit_${DATE}/results

for bam in $allbam
do
    leukid=$(echo $bam | awk 'BEGIN {FS="/"}{print $15}' | awk 'BEGIN {FS=".bam"}{print $1}')
    bsub -We 59 -n $CORES -R 'rusage[mem=16]' -J cnvkit_$leukid -eo cnvkit_err_$leukid -oo cnvkit_out_$leukid ./cnvkit_onesample.sh $DATE $leukid $bam $CORES
done
