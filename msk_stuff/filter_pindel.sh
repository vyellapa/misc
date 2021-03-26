rm run_pindel

PROJ="256"

lk get_outdirs -fi name PINDEL -fi projects__pk ${PROJ} -fi status SUCCEEDED |awk '{print "rsync -aP "$1"/*.flagged.annot.vcf.gz ./"}'|bash

IGH="/ifs/res/leukgen/home/yellapav/record/mm_baits/beds/MM_IGH_intron_nonCHR.bed"
TARGETS="/ifs/res/leukgen/home/yellapav/record/mm_baits/beds/MM_merged_targets.bed"
PINDEL_ANNO="/ifs/res/leukgen/projects/219/ANALYSES/annot_pindel__cosmic_81_exac_03__12/annot_pindel__cosmic_81_exac_03__12.tgd.unfiltered.tsv.gz"


mkdir -p logs

for i in `ls *I-H-100540-N1-1-D1-2.*.vcf.gz`
do

SNAME=`echo ${i} |awk -F'_vs_' '{print $1}'`
NNAME=`echo ${i} |awk -F'_vs_|.fla' '{print $2}'`

echo ${SNAME}

less ${i} |java -jar /opt/common/CentOS_6/snpEff/snpEff_v4.1a/SnpSift.jar intervals ${TARGETS}  | java -jar /opt/common/CentOS_6/snpEff/snpEff_v4.1a/SnpSift.jar filter "(GEN[1].PP>1) & (GEN[1].NP>1) & (FILTER = 'PASS')" > temp
less ${i} |java -jar /opt/common/CentOS_6/snpEff/snpEff_v4.1a/SnpSift.jar intervals ${TARGETS}  | java -jar /opt/common/CentOS_6/snpEff/snpEff_v4.1a/SnpSift.jar filter "(GEN[1].PP>1) & (GEN[1].NP>1) " > temp
cat temp |java -jar /opt/common/CentOS_6/snpEff/snpEff_v4.1a/SnpSift.jar intervals ${IGH} > ${SNAME}.pindel.igh.vcf 
cat temp |java -jar /opt/common/CentOS_6/snpEff/snpEff_v4.1a/SnpSift.jar intervals ${IGH} -x > ${SNAME}.pindel.filter.vcf

mkdir -p MAF
mkdir -p VCF
echo "bsub -oo logs/${SNAME}.log cmo_vcf2maf --ncbi-build GRCh37 --output-maf ${SNAME}.pindel.filter.maf --input-vcf ${SNAME}.pindel.filter.vcf --tumor-id ${SNAME} --normal-id ${NNAME} --vep-release 86 --version 1.6.11" >> run_pindel
echo "bsub -oo logs/${SNAME}.igh.log cmo_vcf2maf --ncbi-build GRCh37 --output-maf ${SNAME}.pindel.igh.maf --input-vcf ${SNAME}.pindel.igh.vcf --tumor-id ${SNAME} --normal-id ${NNAME} --vep-release 86 --version 1.6.11" >> run_pindel

done



cat *.pindel.filter.vcf > temp

mv -f *.maf ./MAF/
mv -f *.vcf ./VCF/
mv -f *.flagged.annot.vcf.gz ./VCF/

#Annotate with MMRF and normal frequency
python ~/local/msk_scripts/annotate_mmrf_v4.py annotate_tsv_freq --in_tsv_gz /ifs/res/leukgen/projects/${PROJ}/ANALYSES/annot_pindel__cosmic_81_exac_03__12/annot_pindel__cosmic_81_exac_03__12.tgd.annot.bed.tsv.gz --annotation_tsv ${PINDEL_ANNO} > p${PROJ}_pindel_mmrf.normals.tsv


## Remove IGH calls
perl -e '%s=map{chomp;@a=split(/\t/);($a[2],1)}`less $ARGV[0]`; foreach $l (`less $ARGV[1]`) {chomp($l); @a=split(/\t/,$l);if($l=~m/^ID_VARIANT/) {print "$l\n"; next;} if($s{$a[0]}) {print "$l\n";} }' temp p${PROJ}_pindel_mmrf.normals.tsv > p${PROJ}_P_nonIGH

less p${PROJ}_P_nonIGH | awk -F'\t' '{if($22!="synonymous_codon" && $22!~/stream_variant/ && $22!="NA" && $22!~/UTR/ && $22!~/intron/ && $22!="nc_transcript_variant") print $0}' > p${PROJ}_P_nonIGH_NS


#Germline SNP filter
less p${PROJ}_P_nonIGH_NS | awk -F'\t' '{OFS="\t"; flag=0;for(i=56;i<63;i++) {if($i>0.03 && $i!="NA") {flag=1;}} if($70>0.03 && $70!="NA") {flag=1;} if(flag==0 || NR==1 ) {print $0}}' > p${PROJ}_P_SNP3

less p${PROJ}_P_SNP3 | awk -F'\t' '{OFS="\t"; flag=0;for(i=56;i<63;i++) {if($i>0.001 && $i!="NA" && $23!="PASS") {flag=1;}} if($70>0.001 && $70!="NA" && $23!="PASS") {flag=1;} if(flag==0 || NR==1 || $55~/^GENOMIC_EXACT/ || $55~/^GENOMIC_POS/) {print $0}}' > p${PROJ}_P_SNP01

#not pass and not genomic
less p${PROJ}_P_SNP01 | awk -F'\t' '{OFS="\t"; flag=0; if($23!="PASS" && $71=="NA" && $77=="NA") {flag=1;} if(flag==0 || NR==1 || $55~/^GENOMIC_EXACT/ || $55~/^GENOMIC/) {print $0}}' > p${PROJ}_P_MMRF_BOL

#If found in atleast 4 nomals
less p${PROJ}_P_MMRF_BOL | awk -F'\t' '{OFS="\t"; flag=0; if($81!="NA" && $81>=1) {flag=1;} if(flag==0 || NR==1) {print $0} }' > p${PROJ}_P_NORMAL1

#not PASS or genomic or MMRF and is found in 1 normal
less p${PROJ}_P_NORMAL1 | awk -F'\t' '{OFS="\t"; flag=0; if($23!="PASS" && $71!="NA" && $77!="NA" && $81!="NA" && $81>0) {flag=1;} if(flag==0 || NR==1 || $55~/^GENOMIC_EXACT/ || $55~/^GENOMIC_POS/) {print $0} }' > p${PROJ}_P_UPLOAD



