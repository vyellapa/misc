echo "# chr1,start1,end1,chr2,start2,end2,id/name,brass_score,strand1,strand2,repeats,np_sample_count,tumour_count,normal_count,np_count,bkpt_distance,sample,sample_type,names,count,bal_trans,inv,occL,occH,copynumber_flag,range_blat" |tr ',' '\t' > t

for i in `ls I-H-1*/intermediates/*groups.gz`
do

NAME=`echo ${i} |awk -F'/' '{print $NF}' | awk -F'_50_vs_' '{print $1}'`
less ${i} | grep "^#" > ${NAME}.groups
less t | grep "^#" >> ${NAME}.groups
less ${i} |grep -v "^#"| awk -v n=${NAME} -F'\t' '{OFS="\t";if($10>5) {print $1,$3,$4,$5,$7,$8,$9,$10, $2, $6,"0","0", "0","0","0","-1",n,"T", $113, "67\t1058\t0\t0\t2\t0\t0"}}' >> ${NAME}.groups

/ifs/work/leukgen/bin/perlbrew/0.74/perls/5.18.4/bin/perl /ifs/work/leukgen/opt/cgp/5.18.4/grass/1.1.6/bin/grass.pl -genome_cache /ifs/work/leukgen/ref/homo_sapiens/GRCh37d5/vagrent/Homo_sapiens_KnC.GRCh37.75.vagrent.cache.gz -ref /ifs/work/leukgen/ref/homo_sapiens/GRCh37d5/genome/gr37.fasta -species HUMAN -assembly GRCH37D5 -platform CARLOS -protocol WGS -tumour I-H-121266-T1-3-D1-1 -normal I-H-121266-N1-2-D1-1 -file ${NAME}.groups

done


