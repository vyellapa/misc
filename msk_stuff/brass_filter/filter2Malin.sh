PROJ="256"
TSTAMP=`date|awk '{OFS="_";print $2,$3,$NF}'`
perl ~/brass_repeatMask.pl E-H-133050-T1-1-D1-1.filtered |head -n1 | cut -f1,2,4,5,8,13,14,16,17,20,27-30,32,36-39,41,45,48,49,50 > p${PROJ}_brassR3_50bp_results_${TSTAMP}.tsv

for i in `ls *.filtered`
do


perl ~/brass_repeatMask.pl ${i} > temp
cat temp | cut -f1,2,4,5,8,13,14,16,17,20,27-30,32,36-39,41,45,48,49,50 >> p${PROJ}_brassR3_50bp_results_${TSTAMP}.tsv


done


