#Usage: ./elli3_abPairs.sh a2.bedpe /ifs/res/leukgen/local/opt/leukdc/data/workflows/25/75/12575/data/bam/E-H-109099-T2-1-D1-1.bam > T2_trans_abPairs.bedpe 
PROJ="162"
for z in `ls p162_delly_results_sep28.tsv` 
do


cat ${z} | grep CHR | awk '{print $0"\tMAPQ1_50\tMATCH1_50\tSUPP1_50\tMAPQ2_50\tMATCH2_50\tSUPP2_50"}' > alll_50


#Usage: ./elli3_abPairs.sh a2.bedpe /ifs/res/leukgen/local/opt/leukdc/data/workflows/25/75/12575/data/bam/E-H-109099-T2-1-D1-1.bam > T2_trans_abPairs.bedpe 

for i in `cat ${z} | tr '\t' '!'|sed 's/ //g'|awk -F'!' '{if($7=="TRA") {print $0}}'`
do



GR=`echo ${i} |awk -F'!' '{print $13}'`

BAM=`cat bam_paths_50 |grep ${GR} |awk -F'\t' '{print $1}'`
echo "Processing translocations ${GR} ..."
CHR1=`echo ${i} |awk -F'!' '{print $1}'`
CHR2=`echo ${i} |awk -F'!' '{print $4}'`

ST1=`echo ${i} |awk -F'!' '{print $2-1000}'`
SP1=`echo ${i} |awk -F'!' '{print $2+1000}'`

ST2=`echo ${i} |awk -F'!' '{print $5-1000}'`
SP2=`echo ${i} |awk -F'!' '{print $5+1000}'`



echo ${i} > m1
echo ${i} > m2
#echo "samtools view ${BAM} ${CHR1}:${ST1}-${SP1} | awk -v st2=${ST2} -v chr1=${CHR1} -v chr2=${CHR2} -F'\t' '{if("\$3"==chr1 && "\$7"==chr2 && (st2-4000)<"\$8" && (st2+5000)>"\$8") {print "\$0"}}'"
samtools view ${BAM} ${CHR1}:${ST1}-${SP1} | awk -v st2=${ST2} -v chr1=${CHR1} -v chr2=${CHR2} -F'\t' '{if($3==chr1 && $7==chr2 && (st2-4000)<$8 && (st2+5000)>$8) {print $0}}' >> m1
samtools view ${BAM} ${CHR2}:${ST2}-${SP2} | awk -v st1=${ST1} -v chr1=${CHR2} -v chr2=${CHR1} -F'\t' '{if($3==chr1 && $7==chr2 && (st1-4000)<$8 && (st1+5000)>$8) {print $0}}' >> m2




perl -e '$match=0;$sum=0;$row=-1;foreach $l (`cat m1`) {++$row;next if($row eq 0); chomp($l); @a=split(/\t/,$l);$sum+=$a[4];@m=$a[5]=~m/(\d+)M/; map{$match += $_ }@m;} if($row>0) {$sum=$sum/$row;$match=$match/$row;print "$sum\t$match\t$row\n";} else {print "0\t0\t0";}' > temp
MAPQ1=`cat temp|cut -f1`
MATCH1=`cat temp|cut -f2`
SUPP1=`cat temp|cut -f3`


perl -e '$match=0;$sum=0;$row=-1;foreach $l (`cat m2`) {++$row;next if($row eq 0); chomp($l); @a=split(/\t/,$l);$sum+=$a[4];@m=$a[5]=~m/(\d+)M/; map{$match += $_ }@m;} if($row>0) {$sum=$sum/$row;$match=$match/$row; print "$sum\t$match\t$row\n";} else {print "0\t0\t0";}' > temp

MAPQ2=`cat temp|cut -f1`
MATCH2=`cat temp|cut -f2`
SUPP2=`cat temp|cut -f3`



echo "${i}	${MAPQ1}	${MATCH1}	${SUPP1}	${MAPQ2}	${MATCH2}	${SUPP2}" |tr '!' '\t' >> alll_50


done



#cat ${NAME} |awk -F'\t' '{if(($47>10 && $50>10 && $48>50 && $51>50 && $49<800 && $54<800) || $47=="NA" || $47=="MAPQ1") {print $0}}' > ${NAME1}

rm -f m1 m2 temp


done