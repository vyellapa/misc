for i in `cat ../delly/bam_paths`
do

NAME=`echo ${i} |awk -F'/' '{print $NF}' | awk -F'.' '{print $1}'`

#ln -s ${i} ./${NAME}.bam
#ln -s ${i}.bai ./${NAME}.bam.bai
#java -jar /home/yellapav/local/bin/gridss-2.0.0-gridss-jar-with-dependencies.jar ASSEMBLY=I-H-135327-T1-1-D1-1.assembly.bam R=/ifs/work/leukgen/ref/homo_sapiens/GRCh37d5/genome/gr37.fasta I=I-H-135327-T1-1-D1-1.bam O=I-H-135327-T1-1-D1-1_output BL=/home/yellapav/local/bin/girids_exclude.bed


#ava -jar /home/yellapav/local/bin/gridss-2.0.0-gridss-jar-with-dependencies.jar ASSEMBLY=I-H-135327-T1-1-D1-1.assembly.bam R=/ifs/work/leukgen/ref/homo_sapiens/GRCh37d5/genome/gr37.fasta I=I-H-135327-T1-1-D1-1.bam O=I-H-135327-T1-1-D1-1_output.vcf BL=/home/yellapav/local/bin/girids_exclude.bed
echo "bsub -oo logs/${NAME}.log -We 1:59 -n 4 -R "rusage[mem=20]" java -jar /home/yellapav/local/bin/gridss-2.0.0-gridss-jar-with-dependencies.jar ASSEMBLY=${NAME}.assembly.bam R=/ifs/work/leukgen/ref/homo_sapiens/GRCh37d5/genome/gr37.fasta I=${NAME}.bam O=${NAME}_output.vcf BL=/home/yellapav/local/bin/girids_exclude.bed THREADS=4"

done
