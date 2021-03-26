#!/bin/bash

DATE=$(date -d '1 hour ago' '+%Y-%m-%d')

CORES=8

BAIT_BED=/ifs/res/papaemme/users/gg10/iwg-pm/bait-design/iwg-pm-targets.bed
ACCESS_BED=/ifs/work/leukgen/ref/homo_sapiens/GRCh37d5/gr37.access.bed
FASTA=/ifs/work/leukgen/ref/homo_sapiens/GRCh37d5/gr37.fasta
OUT=cnvkit_${DATE}
OUT_REF=${OUT}/reference.cnn
OUT_DIR=${OUT}/results

mkdir -p $OUT


#1) symlink the Normal bam files
echo "created symlinks..."
mybam=$(lk get_data -fi projects 105 -fi specimen__source_type NORMAL)

for bam in $mybam
do
	leukid=$(echo $bam | awk 'BEGIN {FS="bam/"}{print $2}' | awk 'BEGIN {FS=".bam"}{print $1}')
	echo $leukid
	ln -s $bam ${leukid}.Normal.bam
done

echo "symlinks created!"

#3) take a random tumor sample

ln -s "/ifs/res/leukgen/local/opt/leukdc/data/workflows/66/98/26698/data/bam/E-H-121171-T1-1-D1-1.bam" tmp.Tumor.bam

#2) run cnvkit for one tumor sample and build the reference 

echo "lauching cnvkit..."

cnvkit.py batch *Tumor.bam --normal *Normal.bam --targets $BAIT_BED --access $ACCESS_BED --fasta $FASTA --output-reference $OUT_REF --output-dir $OUT_DIR -p $CORES

#3) clean output

rm *Normal.bam* tmp.Tumor.bam* $OUT/results/tmp*
