ONE=`ls *.filtered|head -n1`
perl ~/brass_repeatMask.pl ${ONE} |head -n1 | cut -f1,2,4,5,8,13,14,16,17,20,27-30,32,36-39,41,45,48,49,50 > p251_brassR3_50bp_results_oct19.tsv

for i in `ls *.filtered`
do


perl ~/brass_repeatMask.pl ${i} > temp
cat temp | cut -f1,2,4,5,8,13,14,16,17,20,27-30,32,36-39,41,45,48,49,50 |grep -v "^chr1" >> p251_brassR3_50bp_results_oct19.tsv


done


