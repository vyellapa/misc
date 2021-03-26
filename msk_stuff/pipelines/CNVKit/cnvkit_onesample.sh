#!/bin/bash
set -e
set -u
set -o pipefail
#echo $HOSTNAME
#echo $LSB_JOBID

# DATE=$(date -d '1 hour ago' '+%Y-%m-%d')
# DATE=2016-10-27

DATE=$1
LEUKID=$2
BAM=$3
CORES=$4

CUR_DIR=`pwd`

OUT=cnvkit_${DATE}
echo "out dir is $OUT"
OUT_REF=${OUT}/reference.cnn
OUT_DIR=${OUT}/results
OUT_DIR_TMP=${OUT}/results_$LEUKID

#cd $CUR_DIR

mkdir -p $OUT_DIR

#1) symlink the bam files
echo "created symlinks..."
ln -s $BAM $LEUKID.Tumor.bam
echo "symlinks created!"

#2) run cnvkit 
echo "lauching cnvkit..."
cnvkit.py batch $LEUKID.Tumor.bam --reference $OUT_REF --output-dir $OUT_DIR_TMP --diagram --scatter -p $CORES
echo "done!"

#3) clean output
echo "clean"
mv $OUT_DIR_TMP/$LEUKID* $OUT_DIR
rm -r $OUT_DIR_TMP
rm $LEUKID.Tumor.bam
rm $LEUKID.Tumor.bam.bai
echo "finished!"
