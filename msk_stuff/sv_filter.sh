#Usage: ./elli3_abPairs.sh a2.bedpe /ifs/res/leukgen/local/opt/leukdc/data/workflows/25/75/12575/data/bam/E-H-109099-T2-1-D1-1.bam > T2_trans_abPairs.bedpe 
for z in `ls ../../E-H-109*/brass_*/*.annot.bedpe`
do

echo "Sample ${z}"
GR=`echo ${z} |awk -F'/' '{print $NF}' |awk -F'_vs_' '{print $1}'` 


SV=${z}
BAM=`lk get data -p 163 -v |grep ${GR} |awk -F'\t' '{print $2}'`
NAME=`echo ${SV} | awk -F'/' '{print $NF}' |sed 's/annot/annot.filter/g'`
NAME1=`echo ${SV} | awk -F'/' '{print $NF}' |sed 's/annot/annot.filtered/g'`

cat ${SV} | grep chr  |awk '{print $0"\tMAPQ1\tMATCH1\tSUPP1\tMAPQ2\tMATCH2\tSUPP2"}' > ${NAME}

for i in `cat ${SV} | tr '\t' '!'|sed 's/ //g'|awk -F'!' '{if($12=="deletion") {print $0}}'`
do

echo "Processing deletions .."
CHR1=`echo ${i} |awk -F'!' '{print $1}'`
CHR2=`echo ${i} |awk -F'!' '{print $4}'`

ST1=`echo ${i} |awk -F'!' '{print $2-1000}'`
SP1=`echo ${i} |awk -F'!' '{print $2+1000}'`

ST2=`echo ${i} |awk -F'!' '{print $5-1000}'`
SP2=`echo ${i} |awk -F'!' '{print $5+1000}'`



echo ${i} > m1
echo ${i} > m2
samtools view ${BAM} ${CHR1}:${ST1}-${SP1} | awk -v st2=${ST2} -v chr1=${CHR1} -v chr2=${CHR2} -F'\t' '{if($3==chr1 && $7=="=" && st2<$8 && (st2+2000)>$8) {print $0}}' >> m1
samtools view ${BAM} ${CHR2}:${ST2}-${SP2} | awk -v st1=${ST1} -v chr1=${CHR2} -v chr2=${CHR1} -F'\t' '{if($3==chr1 && $7=="=" && (st1)<$8 && (st1+2000)>$8) {print $0}}' >> m2




perl -e '$match=0;$sum=0;$row=-1;foreach $l (`cat m1`) {++$row;next if($row eq 0); chomp($l); @a=split(/\t/,$l);$sum+=$a[4];@m=$a[5]=~m/(\d+)M/; map{$match += $_ }@m;} if($row>0) {$sum=$sum/$row;$match=$match/$row;print "$sum\t$match\t$row\n";} else {print "0\t0\t0";}' > temp
MAPQ1=`cat temp|cut -f1`
MATCH1=`cat temp|cut -f2`
SUPP1=`cat temp|cut -f3`


perl -e '$match=0;$sum=0;$row=-1;foreach $l (`cat m2`) {++$row;next if($row eq 0); chomp($l); @a=split(/\t/,$l);$sum+=$a[4];@m=$a[5]=~m/(\d+)M/; map{$match += $_ }@m;} if($row>0) {$sum=$sum/$row;$match=$match/$row; print "$sum\t$match\t$row\n";} else {print "0\t0\t0";}' > temp

MAPQ2=`cat temp|cut -f1`
MATCH2=`cat temp|cut -f2`
SUPP2=`cat temp|cut -f3`


echo "${i}	${MAPQ1}	${MATCH1}	${SUPP1}	${MAPQ2}	${MATCH2}	${SUPP2}" |tr '!' '\t' >> ${NAME}




done
#Usage: ./elli3_abPairs.sh a2.bedpe /ifs/res/leukgen/local/opt/leukdc/data/workflows/25/75/12575/data/bam/E-H-109099-T2-1-D1-1.bam > T2_trans_abPairs.bedpe 

for i in `cat ${SV} | tr '\t' '!'|sed 's/ //g'|awk -F'!' '{if($12=="translocation") {print $0}}'`
do

echo "Processing translocations ..."
CHR1=`echo ${i} |awk -F'!' '{print $1}'`
CHR2=`echo ${i} |awk -F'!' '{print $4}'`

ST1=`echo ${i} |awk -F'!' '{print $2-1000}'`
SP1=`echo ${i} |awk -F'!' '{print $2+1000}'`

ST2=`echo ${i} |awk -F'!' '{print $5-1000}'`
SP2=`echo ${i} |awk -F'!' '{print $5+1000}'`



echo ${i} > m1
echo ${i} > m2
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



echo "${i}	${MAPQ1}	${MATCH1}	${SUPP1}	${MAPQ2}	${MATCH2}	${SUPP2}" |tr '!' '\t' >> ${NAME}


done


cat ${SV} | awk -F'\t' '{if($12!="translocation" && $12!="deletion" && $1!~/^#/) {print $0"\tNA\tNA\tNA\tNA\tNA\tNA";}}' >> ${NAME}

cat ${NAME} |awk -F'\t' '{if(($47>10 && $50>10 && $48>50 && $51>50 && $49<800 && $54<800) || $47=="NA" || $47=="MAPQ1") {print $0}}' > ${NAME1}

rm -f m1 m2 temp


done
