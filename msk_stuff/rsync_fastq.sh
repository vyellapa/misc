lk get_data --dtype fastq -fi project 292 -v > fastqs_all

mkdir -p logs
for i in `cat fastqs_all |awk -F'\t' '{print $1"!"$2 }'`
do
NAME=`echo ${i} | awk -F'!' '{print $1}'`
FQ=`echo ${i} | awk -F'!' '{print $2}'`
FQ_PATH=`readlink -f ${FQ}`
FQ_NAME=`echo ${FQ} | awk -v n=${NAME} -F'/' '{print n"_"$NF}'`
FQ_NAME_50=`echo ${FQ} | awk -v n=${NAME} -F'/' '{print n"_50bp_"$NF}'`

## Oldway Below ##
#echo "rsync -aP ${FQ_PATH} fastq_100/${FQ_NAME}"
#echo "bsub -oo logs/${FQ_NAME_50}.log \"less fastq_100/${FQ_NAME} | /home/yellapav/local/downloads/fastx/bin/fastx_trimmer -t 50 -z -Q33 -o ${FQ_NAME_50}\""

echo "bsub -oo logs/${FQ_NAME_50}.log \"less ${FQ_PATH} | /home/yellapav/local/downloads/fastx/bin/fastx_trimmer -t 51 -z -Q33 -o ${FQ_NAME_50}\""



done
