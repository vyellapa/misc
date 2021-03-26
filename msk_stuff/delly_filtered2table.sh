PROJ="220"

echo "1" |awk -F'\t' '{OFS="\t"; print "CHR1","START1","END1","CHR2", "START2", "END2", "SVTYPE", "GENE","TUMOR_JUNCTIONS_READS","TUMOR_JUNCTIONS_PAIRS","NORMAL_JUNCTIONS_READS","NORMAL_JUNCTION_PAIRS","SAMPLE" }' > p${PROJ}_delly_results.tsv



for i in `ls *.vcf`
do

NAME=`echo ${i} |awk -F'.' '{print $1}'`

cat ${i} |sed 's/ID=END/ID=END2/g'|sed 's/;END=/;END2=/g' > temp.vcf
java -jar /opt/common/CentOS_6/snpEff/snpEff_v4.1a/SnpSift.jar extractFields temp.vcf CHROM POS CHR2 END2 SVTYPE SVMETHOD GEN[0].DV GEN[0].RV GEN[1].DV GEN[1].RV| sed 's/EMBL.DELLYv0.7.6,//g'|sed 's/EMBL.DELLYv0.7.6/NA/g' |sed 's/GEN\[0\]\.DV/TUMOR_JUNCTIONS_PAIRS/g'|sed 's/GEN\[0\]\.RV/TUMOR_JUNCTIONS_READS/g' |sed 's/GEN\[1\]\.DV/NORMAL_JUNCTIONS_PAIRS/g' |sed 's/GEN\[1\]\.RV/NORMAL_JUNCTION_READS/g' > temp

cat temp |grep -v "^#" | awk -v n=${NAME} -F'\t' '{OFS="\t"; print $1, $2, $2+100, $3, $4,$4+100, $5,$6,$7,$8,$9,$10,n}'  |awk -F'\t' '{flag=0; if($7=="TRA" && $1!=14 && $4!=14) {flag=1;} if(flag==0) {print $0} }' >> p${PROJ}_delly_results.tsv



done
