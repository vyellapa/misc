for i in `ls I-H*.gz| awk -F'/' '{print $NF}' | awk -F'_50bp_' '{print $1}'|sort -u`
do

#NAME=`echo ${i} | awk -F'/' '{print $NF}' | awk -F'_50bp_' '{print $1}'`
NAME=${i}

echo "mkdir -p ${NAME}"
echo "mv ${NAME}*.gz ./${NAME}"

echo "cd ${NAME}"
echo "bsub -We 1:59 -oo bwa.log -n 4 -R \"rusage[mem=12]\" bwa_mem.pl -t 4 -r /home/yellapav/local/resources/hs37d5.fa -o myout -s ${NAME} ${NAME}*_[12].fastq.gz"
echo "cd -"


done
