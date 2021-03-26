GENOME="/ifs/work/leukgen/ref/homo_sapiens/GRCh37d5/genome/gr37.fasta"
GTF="/home/yellapav/local/resources/Homo_sapiens.GRCh37.75.gtf"
GDIR="/home/yellapav/local/resources/star_h37d5_ens75"
TEMP="/ifs/work/leukgen/home/yellapav/star_align/tmp"

CDS="/home/yellapav/local/resources/Homo_sapiens.GRCh37.75.CDS.bed"
RIBO="/home/yellapav/local/resources/RIBOSOMAL.bed"
CDS="/home/yellapav/local/resources/refFlat.intervals" 
RIBO="/home/yellapav/local/resources/RIBOSOMAL.intervals"

DB_129="/ifs/res/leukgen/home/yellapav/resources/gatk_bundle/ftp.broadinstitute.org/bundle/2.8/b37/dbsnp_138.b37.excluding_sites_after_129.vcf.gz"
MY_CDS="/ifs/res/leukgen/home/yellapav/record/mm_baits/myType/beds/cds_merged.interval_list"
CYTO="/ifs/res/leukgen/home/yellapav/record/mm_baits/myType/beds/cyto_beds.interval_list"
FINGER="/ifs/res/leukgen/home/yellapav/record/mm_baits/myType/beds/fingerPrint.interval_list"
FOCAL="/ifs/res/leukgen/home/yellapav/record/mm_baits/myType/beds/focal.interval_list"
IGH="/ifs/res/leukgen/home/yellapav/record/mm_baits/myType/beds/MM_IGH_intron.interval_list"
MY_ALL="/ifs/res/leukgen/home/yellapav/record/mm_baits/myType/beds/MM_MSK_permiss_capture_targets.interval_list"




BAM=$1
NAME=`echo ${BAM} |awk -F'/' '{print $NF}' | awk -F'.bam' '{print $1}'`

mkdir -p logs

bsub -oo logs/${NAME}_CollectAlignmentSummaryMetrics.log -R "rusage[mem=8]" "java -jar /opt/common/CentOS_6-dev/picard/v1.140/picard.jar CollectAlignmentSummaryMetrics INPUT=${BAM} OUTPUT=${NAME}.picard.alignment_summary_metrics.tsv REFERENCE_SEQUENCE=${GENOME} VALIDATION_STRINGENCY=LENIENT"

bsub -oo logs/${NAME}_CollectInsertMetrics.log -R "rusage[mem=8]" "java -jar /opt/common/CentOS_6-dev/picard/v1.140/picard.jar CollectInsertSizeMetrics INPUT=${BAM} OUTPUT=${NAME}.picard.insertSize_metrics.tsv H=${NAME}.insert.size.histogram.pdf VALIDATION_STRINGENCY=LENIENT"


## For WGS only
#bsub -oo logs/${NAME}_CollectRawWgsMetrics.log -R "rusage[mem=48]" "java -jar /opt/common/CentOS_6-dev/picard/v1.140/picard.jar CollectRawWgsMetrics INPUT=${BAM} OUTPUT=${NAME}.picard.raw_wgs_metrics.tsv REFERENCE_SEQUENCE=${GENOME} INCLUDE_BQ_HISTOGRAM=true VALIDATION_STRINGENCY=LENIENT"

#bsub -oo logs/${NAME}_CollectWgsMetrics.log -R "rusage[mem=48]" "java -jar /opt/common/CentOS_6-dev/picard/v1.140/picard.jar CollectWgsMetrics INPUT=${BAM} OUTPUT=${NAME}.picard.wgs_metrics.tsv REFERENCE_SEQUENCE=${GENOME} VALIDATION_STRINGENCY=LENIENT"

#bsub -oo logs/${NAME}_CollectWgsMetricsNonZero.log -R "rusage[mem=48]" "java -jar /opt/common/CentOS_6-dev/picard/v1.140/picard.jar CollectWgsMetricsWithNonZeroCoverage INPUT=${BAM} OUTPUT=${NAME}.picard.nonZero_wgs_metrics.tsv CHART=${NAME}.picard.nonZero.pdf REFERENCE_SEQUENCE=${GENOME} VALIDATION_STRINGENCY=LENIENT"

####################


bsub -oo logs/${NAME}_CollectSeqArticats.log -R "rusage[mem=8]" "java -jar /opt/common/CentOS_6-dev/picard/v1.140/picard.jar CollectSequencingArtifactMetrics INPUT=${BAM} OUTPUT=${NAME}.picard.artifact_metrics.tsv REFERENCE_SEQUENCE=${GENOME} VALIDATION_STRINGENCY=LENIENT"


bsub -oo logs/${NAME}_overall_hsMetrics.log -R "rusage[mem=8]" "java -jar /opt/common/CentOS_6-dev/picard/v1.140/picard.jar CalculateHsMetrics INPUT=${BAM} OUTPUT=${NAME}.overall_hs_metrics.tsv BAIT_INTERVALS=${MY_ALL} TARGET_INTERVALS=${MY_ALL} REFERENCE_SEQUENCE=${GENOME} VALIDATION_STRINGENCY=LENIENT"
bsub -oo logs/${NAME}_cds_hsMetrics.log -R "rusage[mem=8]" "java -jar /opt/common/CentOS_6-dev/picard/v1.140/picard.jar CalculateHsMetrics INPUT=${BAM} OUTPUT=${NAME}.cds_hs_metrics.tsv BAIT_INTERVALS=${MY_CDS} TARGET_INTERVALS=${MY_CDS} REFERENCE_SEQUENCE=${GENOME} VALIDATION_STRINGENCY=LENIENT"
bsub -oo logs/${NAME}_cyto_hsMetrics.log -R "rusage[mem=8]" "java -jar /opt/common/CentOS_6-dev/picard/v1.140/picard.jar CalculateHsMetrics INPUT=${BAM} OUTPUT=${NAME}.cyto_hs_metrics.tsv BAIT_INTERVALS=${CYTO} TARGET_INTERVALS=${CYTO} REFERENCE_SEQUENCE=${GENOME} VALIDATION_STRINGENCY=LENIENT"
bsub -oo logs/${NAME}_finger_hsMetrics.log -R "rusage[mem=8]" "java -jar /opt/common/CentOS_6-dev/picard/v1.140/picard.jar CalculateHsMetrics INPUT=${BAM} OUTPUT=${NAME}.finger_hs_metrics.tsv BAIT_INTERVALS=${FINGER} TARGET_INTERVALS=${FINGER} REFERENCE_SEQUENCE=${GENOME} VALIDATION_STRINGENCY=LENIENT"
bsub -oo logs/${NAME}_focal_hsMetrics.log -R "rusage[mem=8]" "java -jar /opt/common/CentOS_6-dev/picard/v1.140/picard.jar CalculateHsMetrics INPUT=${BAM} OUTPUT=${NAME}.focal_hs_metrics.tsv BAIT_INTERVALS=${FOCAL} TARGET_INTERVALS=${FOCAL} REFERENCE_SEQUENCE=${GENOME} VALIDATION_STRINGENCY=LENIENT"
bsub -oo logs/${NAME}_igh_hsMetrics.log -R "rusage[mem=8]" "java -jar /opt/common/CentOS_6-dev/picard/v1.140/picard.jar CalculateHsMetrics INPUT=${BAM} OUTPUT=${NAME}.igh_hs_metrics.tsv BAIT_INTERVALS=${IGH} TARGET_INTERVALS=${IGH} REFERENCE_SEQUENCE=${GENOME} VALIDATION_STRINGENCY=LENIENT"


bsub -oo logs/${NAME}_CollectOxoGMetrics.log -R "rusage[mem=8]" "java -jar /opt/common/CentOS_6-dev/picard/v1.140/picard.jar CollectOxoGMetrics INPUT=${BAM} OUTPUT=${NAME}.CollectOxoGMetrics.tsv REFERENCE_SEQUENCE=${GENOME} DB_SNP=${DB_129} VALIDATION_STRINGENCY=LENIENT"


bsub -oo logs/${NAME}_MarkDuplicates.log -R "rusage[mem=8]" "java -Xmx4G -jar /opt/common/CentOS_6-dev/picard/v1.140/picard.jar MarkDuplicates INPUT=${BAM} OUTPUT=${NAME}.mdup.bam M=${NAME}.mdup.txt VALIDATION_STRINGENCY=LENIENT"


bsub -oo logs/${NAME}_samstats.log "samtools flagstat ${BAM} > ${NAME}.samStats"



