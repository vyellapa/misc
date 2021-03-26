mkdir -p logs

for i in `less bam_paths`
do

NAME=`echo ${i} | awk -F'/' '{print $NF}' |awk -F'.' '{print $1}'`

N_BAM="/ifs/res/leukgen/projects/222/RESULTS/brass_50/bams_normal/I-H-100540-N1-1-D1-2/myout/I-H-100540-N1-1-D1-2.bam"
SHORT_NAME=`echo ${NAME} | awk '{a=substr($1,1,10); print a}'`
T_BAM=`cat bam_paths | grep ${NAME}`
#N_BAM=`cat bam_paths | grep ${SHORT_NAME}|grep "N1"|head -1`

mkdir -p ${NAME}_brass
echo "bsub -oo logs/${NAME}_brass.log lkcgp brass --tumour-bam ${T_BAM} --normal-bam ${N_BAM} --outdir /ifs/res/leukgen/projects/162/RESULTS/SV/brass/brass/${NAME}_brass"


done
