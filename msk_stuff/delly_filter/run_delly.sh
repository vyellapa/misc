mkdir -p logs
REF="/home/yellapav/local/resources/hs37d5.fa"
EXCL="/home/yellapav/local/downloads/delly/excludeTemplates/human.hg19.excl.tsv"


for i in `cat bam_paths`
do

NBAM="/ifs/res/leukgen/local/opt/leukdc/data/workflows/39/24/23924/data/bam/I-H-100540-N1-1-D1-2.bam"

T1=`echo ${i} |awk -F'/' '{print $NF}' | awk -F'.' '{print $1}'` 
T1_BAM=${i}


echo "bsub -oo logs/${T1}_DEL.log -n 8 -R \"rusage[mem=12]\" /home/yellapav/local/downloads/delly/src/delly call -t DEL -q 10 -s 12 -x ${EXCL} -o ${T1}_DEL.bcf -g ${REF} ${T1_BAM} ${NBAM}"
echo "bsub -oo logs/${T1}_DUP.log -n 8 -R \"rusage[mem=12]\" /home/yellapav/local/downloads/delly/src/delly call -t DUP -q 10 -s 12 -x ${EXCL} -o ${T1}_DUP.bcf -g ${REF} ${T1_BAM} ${NBAM}"
echo "bsub -oo logs/${T1}_INV.log -n 8 -R \"rusage[mem=12]\" /home/yellapav/local/downloads/delly/src/delly call -t INV -q 10 -s 12 -x ${EXCL} -o ${T1}_INV.bcf -g ${REF} ${T1_BAM} ${NBAM}"
echo "bsub -oo logs/${T1}_INS.log -n 8 -R \"rusage[mem=12]\" /home/yellapav/local/downloads/delly/src/delly call -t INS -q 10 -s 12 -x ${EXCL} -o ${T1}_INS.bcf -g ${REF} ${T1_BAM} ${NBAM}"
echo "bsub -oo logs/${T1}_TRA.log -n 8 -R \"rusage[mem=12]\" /home/yellapav/local/downloads/delly/src/delly call -t TRA -q 10 -s 12 -x ${EXCL} -o ${T1}_TRA.bcf -g ${REF} ${T1_BAM} ${NBAM}"
done



