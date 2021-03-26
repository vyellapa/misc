for i in `ls *_TRA.bcf`
do
NAME=`echo ${i} |awk -F'_' '{print $1}'`

bcftools view ${i} | awk -F'\t' '{if($7=="PASS" || $1~/^#/) print $0}' |java -jar /opt/common/CentOS_6/snpEff/snpEff_v3.6/SnpSift.jar filter "( GEN[0].DV > 1 ) & ( GEN[0].RV >= 1 ) & ( GEN[0].GT!=GEN[1].GT )" > temp.vcf
bcftools view ${i} | awk -F'\t' '{if($7=="PASS" || $1~/^#/) print $0}' |java -jar /opt/common/CentOS_6/snpEff/snpEff_v3.6/SnpSift.jar filter "( GEN[0].DV > 0 ) & ( GEN[0].RV > 0 ) " > temp.vcf
python ~/local/msk_scripts/annotate_delly.py annotate_sv_vcf --in_vcf temp.vcf --out_vcf ${NAME}.tra.vcf

done


for i in `ls *_DEL.bcf`
do
NAME=`echo ${i} |awk -F'_' '{print $1}'`

bcftools view ${i} | awk -F'\t' '{if($7=="PASS" || $1~/^#/) print $0}' |java -jar /opt/common/CentOS_6/snpEff/snpEff_v3.6/SnpSift.jar filter "( GEN[0].DV > 0 ) & ( GEN[0].RV > 0 ) & ( GEN[0].GT!=GEN[1].GT ) & ( CHROM != '14' )" > temp.vcf
python ~/local/msk_scripts/annotate_delly.py annotate_sv_vcf --in_vcf temp.vcf --out_vcf ${NAME}.del.vcf

done


for i in `ls *_DUP.bcf`
do
NAME=`echo ${i} |awk -F'_' '{print $1}'`

bcftools view ${i} | awk -F'\t' '{if($7=="PASS" || $1~/^#/) print $0}' |java -jar /opt/common/CentOS_6/snpEff/snpEff_v3.6/SnpSift.jar filter "( GEN[0].DV > 0 ) & ( GEN[0].RV > 0 ) & ( GEN[0].GT!=GEN[1].GT ) & ( CHROM != '14' )" > temp.vcf
python ~/local/msk_scripts/annotate_delly.py annotate_sv_vcf --in_vcf temp.vcf --out_vcf ${NAME}.dup.vcf

done

for i in `ls *_INV.bcf`
do
NAME=`echo ${i} |awk -F'_' '{print $1}'`

bcftools view ${i} | awk -F'\t' '{if($7=="PASS" || $1~/^#/) print $0}' |java -jar /opt/common/CentOS_6/snpEff/snpEff_v3.6/SnpSift.jar filter "( GEN[0].DV > 0 ) & ( GEN[0].RV > 0 ) & ( GEN[0].GT!=GEN[1].GT )" > temp.vcf
python ~/local/msk_scripts/annotate_delly.py annotate_sv_vcf --in_vcf temp.vcf --out_vcf ${NAME}.inv.vcf

done



