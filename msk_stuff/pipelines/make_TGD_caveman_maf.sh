DIR=`pwd`
mkdir -p logs

for i in `ls *.vcf.gz`
do

SNAME=`echo ${i} |awk -F'_vs_' '{print $1}'`
NNAME=`echo ${i} |awk -F'_vs_|.cave' '{print $2}'`
RELAPSE=`echo ${i} | sed 's/-T1-1-D1-1/-T2-1-D1-1/g'`

echo ${SNAME}


less ${i} | awk -F'\t' '{if($1~/^#/ || $7=="PASS" || $7=="VUM;MNP" || $7=="MNP;VUM" || $7=="VUM" || $7=="MNP") {print $0} }' > ${SNAME}.caveman.filter.vcf

mkdir -p MAF
mkdir -p VCF
bsub -oo logs/${SNAME}.log cmo_vcf2maf --ncbi-build GRCh37 --output-maf ${SNAME}.caveman.maf --input-vcf ${SNAME}.caveman.filter.vcf --tumor-id ${SNAME} --normal-id ${NNAME} --version 1.6.11

#mv ${SNAME}.caveman.maf ./MAF/
#mv ${SNAME}.caveman.filter.vcf ./VCF/
done
