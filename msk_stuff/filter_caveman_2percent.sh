rm run_caveman

PROJ="256"

lk get_outdirs -fi name CAVEMAN -fi projects__pk ${PROJ} -fi status SUCCEEDED |awk '{print "rsync -aP "$1"/*.caveman.muts.annot.vcf.gz ./"}'|bash

mkdir -p logs
IGH="/ifs/res/leukgen/home/yellapav/record/mm_baits/beds/MM_IGH_intron_nonCHR.bed"
TARGETS="/ifs/res/leukgen/home/yellapav/record/mm_baits/beds/MM_merged_targets.bed"
#for i in `ls *.caveman.muts.annot.vcf.gz`
for i in `ls *I-H-100540-N1-1-D1-2.caveman.muts.annot.vcf.gz`
do

SNAME=`echo ${i} |awk -F'_vs_' '{print $1}'`
NNAME=`echo ${i} |awk -F'_vs_|.cave' '{print $2}'`
RELAPSE=`echo ${i} | sed 's/-T1-1-D1-1/-T2-1-D1-1/g'`

echo ${SNAME}


less ${i} |java -jar /opt/common/CentOS_6/snpEff/snpEff_v4.1a/SnpSift.jar intervals ${TARGETS}  > temp1
cat temp1 |java -jar /opt/common/CentOS_6/snpEff/snpEff_v4.1a/SnpSift.jar intervals ${IGH} > ${SNAME}.caveman.igh.vcf 
cat temp1 |java -jar /opt/common/CentOS_6/snpEff/snpEff_v4.1a/SnpSift.jar intervals ${IGH} -x > ${SNAME}.caveman.filter.vcf

mkdir -p MAF
mkdir -p VCF
echo "bsub -oo logs/${SNAME}.log cmo_vcf2maf --ncbi-build GRCh37 --output-maf ${SNAME}.caveman.filter.maf --input-vcf ${SNAME}.caveman.filter.vcf --tumor-id ${SNAME} --normal-id ${NNAME} --vep-release 86 --version 1.6.11" >> run_caveman
echo "bsub -oo logs/${SNAME}.igh.log cmo_vcf2maf --ncbi-build GRCh37 --output-maf ${SNAME}.caveman.igh.maf --input-vcf ${SNAME}.caveman.igh.vcf --tumor-id ${SNAME} --normal-id ${NNAME} --vep-release 86 --version 1.6.11" >> run_caveman

done


cat *.caveman.filter.vcf > temp
#~/local/msk_scripts/variant_id.sh tttt ${i} > temp 

#./filter_tsv.sh caveman.txt temp  > a
#cat a | awk -F'\t' '{OFS="\t"; if($22!~/stream/ && $22!="NA" && $22!="synonymous_codon" && $22!~/prime/ && $22!~/intron/) {print $0}}' > caveman.txt


mv -f *.maf ./MAF/
mv -f *.vcf ./VCF/
mv -f *.caveman.muts.annot.vcf.gz ./VCF/



python ~/local/msk_scripts/annotate_mmrf_v4.py annotate_tsv_freq --in_tsv_gz /ifs/res/leukgen/projects/${PROJ}/ANALYSES/annot_caveman__cosmic_81_exac_03__11/annot_caveman__cosmic_81_exac_03__11.tgd.annot.bed.tsv.gz --annotation_tsv /ifs/res/leukgen/projects/219/ANALYSES/annot_caveman__cosmic_81_exac_03__11/annot_caveman__cosmic_81_exac_03__11.tgd.unfiltered.tsv.gz > p${PROJ}_caveman_mmrf.normals.tsv


#UNMATCHED only
cat p${PROJ}_caveman_mmrf.normals.tsv | awk -F'\t' '{if($3=="I-H-100540-N1-1-D1-2" || NR==1) {print $0}}' > p${PROJ}_caveman_mmrf.normals_unmatch.tsv

## Remove IGH calls
perl -e '%s=map{chomp;@a=split(/\t/);($a[2],1)}`less $ARGV[0]`; foreach $l (`less $ARGV[1]`) {chomp($l); @a=split(/\t/,$l);if($l=~m/^ID_VARIANT/) {print "$l\n"; next;} if($s{$a[0]}) {print "$l\n";} }' temp p${PROJ}_caveman_mmrf.normals_unmatch.tsv > p${PROJ}_C_nonIGH

less p${PROJ}_C_nonIGH | awk -F'\t' '{if($22!="synonymous_codon" && $22!~/stream_variant/ && $22!="NA" && $22!~/UTR/ && $22!~/intron/) print $0}' > p${PROJ}_C_nonIGH_NS


#Germline SNP filter
less p${PROJ}_C_nonIGH_NS | awk -F'\t' '{OFS="\t"; flag=0;for(i=74;i<81;i++) {if($i>0.03 && $i!="NA") {flag=1;}} if($88>0.03 && $88!="NA") {flag=1;} if(flag==0 || NR==1) {print $0}}' > p${PROJ}_C_SNP3

less p${PROJ}_C_SNP3 | awk -F'\t' '{OFS="\t"; flag=0;for(i=74;i<81;i++) {if($i>0.001 && $i!="NA") {flag=1;}} if($88>0.001 && $88!="NA") {flag=1;} if(flag==0 || NR==1 || $73~/^GENOMIC_EXACT/ || $73~/^GENOMIC_POS/) {print $0}}' > p${PROJ}_C_SNP01

#not pass and not genomic
less p${PROJ}_C_SNP01 | awk -F'\t' '{OFS="\t"; flag=0; if($23!="PASS" && $89=="NA" && $95=="NA" && $22=="non_synonymous_codon" ) {flag=1;} if(flag==0 || NR==1 || $73~/^GENOMIC_EXACT/ || $73~/^GENOMIC_POS/) {print $0}}' > p${PROJ}_C_MMRF_BOL

#If found in atleast 1 normal
less p${PROJ}_C_MMRF_BOL | awk -F'\t' '{OFS="\t"; flag=0; if($99>=1 && $99!="NA") {flag=1;} if(flag==0 || NR==1 ) {print $0} }' > p${PROJ}_C_NORMAL1

#not PASS or genomic or MMRF and is found in 1 normal
less p${PROJ}_C_NORMAL1 | awk -F'\t' '{OFS="\t"; flag=0; if($23!="PASS" && $89!="NA" && $99>0 && $99!="NA") {flag=1;} if(flag==0 || NR==1 || $73~/^GENOMIC_EXACT/ || $73~/^GENOMIC_POS/) {print $0} }' > p${PROJ}_CAVEMAN_UPLOAD


#If not PASSS and only close/exact in <0.25% of the samples
#cat ttt| awk -F'\t' '{max=0; flag=0; split($90,a,":"); for(i=1;i<=length(a);i++) {if(a[i]>max && $90!="NA") {max=a[i];}} if($23!="PASS" && $89!="NA" && max<2) {flag=1;} if(flag==0 || NR==1 || $73~/^GENOMIC_EXACT/ || $73~/^GENOMIC_POS/) {print $0}   }' > tttt


