#lk get_outdirs -fi name CNVKIT -fi projects__pk 219 -fi status SUCCEEDED | awk '{print "rsync -aP "$1"/*.cns ./"}'|bash

echo "ID	chrom	loc.start	loc.end	num.mark	seg.mean"

for i in `ls *.cns`
do
NAME=`echo ${i} |awk -F'.' '{print $1}'`

cat ${i} |grep -v "^chrom"| awk -v n=${NAME} -F'\t' '{OFS="\t"; print n,$1,$2,$3,$7,$5}' 




done
