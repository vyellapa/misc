#change project number in  the line below

mkdir -p logs

#change project number in  the line below
for i in `lk get_data -fi project 219 |grep -v "I-H-100537-N1-1-D1-2"`
do


N1=`echo ${i}|awk -F'/' '{print $NF}' | awk -F'.' '{print $1}'`
T1=`echo ${i}|awk -F'/' '{print $NF}' | awk -F'.' '{print $1}'|sed 's/-N1-1-D1/-T1-1-D1/g'`

DB="/ifs/res/leukgen/home/yellapav/resources/gatk_bundle/ftp.broadinstitute.org/bundle/2.8/b37/dbsnp_138.b37.vcf.gz"


#Give path to your normal here
N1BAM="/ifs/res/leukgen/local/opt/leukdc/data/workflows/39/22/23922/data/bam/I-H-100537-N1-1-D1-2.bam"

echo "${N1}"

#which snp-pileup
#/home/yellapav/local/bin/snp-pileup
echo "bsub -oo logs/${N1}.log \"snp-pileup -g -q15 -Q20 -P100 -r25,0 ${DB} ${N1}.txt ${N1BAM} ${i}\""



done
