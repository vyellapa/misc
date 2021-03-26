echo "# chr1,start1,end1,chr2,start2,end2,id/name,brass_score,strand1,strand2,repeats,np_sample_count,tumour_count,normal_count,np_count,bkpt_distance,sample,sample_type,names,count,bal_trans,inv,occL,occH,copynumber_flag,range_blat" |tr ',' '\t' > t

for i in `ls E-H*/intermediates/E*.r3`
do

NAME=`echo ${i} |awk -F'/' '{print $NF}' | awk -F'_vs_' '{print $1}'`
GZ=`ls E-H*_brass/intermediates/E*.groups.gz|grep ${NAME}` 
less ${GZ} | grep "^#" > ${NAME}.groups
less t | grep "^#" >> ${NAME}.groups
#less ${i} |grep -v "^#"| awk -v n=${NAME} -F'\t' '{OFS="\t";if($10>5) {print $1,$3,$4,$5,$7,$8,$9,$10, $2, $6,"0","0", "0","0","0","-1",n,"T", $113, "67\t1058\t0\t0\t2\t0\t0"}}' >> ${NAME}.groups

less ${i} |grep -v "^#" >> ${NAME}.groups

/ifs/work/leukgen/bin/perlbrew/0.74/perls/5.18.4/bin/perl /ifs/work/leukgen/opt/cgp/5.18.4/grass/1.1.6/bin/grass.pl -genome_cache /ifs/work/leukgen/ref/homo_sapiens/GRCh37d5/vagrent/Homo_sapiens_KnC.GRCh37.75.vagrent.cache.gz -ref /ifs/work/leukgen/ref/homo_sapiens/GRCh37d5/genome/gr37.fasta -species HUMAN -assembly GRCH37D5 -platform CARLOS -protocol WGS -tumour I-H-121266-T1-3-D1-1 -normal I-H-121266-N1-2-D1-1 -file ${NAME}.groups

done


for i in `ls *_ann.groups`
do

NAME=`echo ${i} |awk -F'_' '{print $1}'`

less ${i} | sed 's/# chr/chr/g' |grep -v "^#" | awk -F'\t' '{if($8>3) print $0}'|awk -F'\t' '{flag=0;if($1==14 && $4==14) {flag=1;} if(flag==0) {print $0;} else {flag=0;}}' > temp

cat temp |awk -F'\t' '{if($1==14 || $4==14 || NR==1) {print $0}}' > temp_
 
less temp_ | awk -F'\t' '{OFS="\t";if($1~/^chr/) {$(NF+1)="svtype";} else if($1!=$4) {$(NF+1)="translocation"} else if($9=="+" && $10=="+") {$(NF+1)="deletion"} else if($9=="-" && $10=="-") {$(NF+1)="tandem-duplication"} else {$(NF+1)="inversion"} print $0}' > ${NAME}.filtered 


##less ${i} | grep -v "^#" | awk -F'\t' '{if($8>20) print $0}' > ${NAME}.filtered 


done 

rm t temp temp_
