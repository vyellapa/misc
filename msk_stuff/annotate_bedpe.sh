for i in `ls -F |grep "^E.*/$"`
do

NAME=`echo ${i} |awk -F'/' '{print $1}'`
SV=`ls ${NAME}/brass*/*.annot.bedpe|grep "${NAME}"`
CNV=`ls ${NAME}/brass*/intermediates/*.ngscn.segments.abs_cn.bg.gz|grep "${NAME}"`


echo "${NAME} ${SV} ${CNV}"
perl -e 'foreach $l (`cat $ARGV[0]`) {chomp($l); @a=split(/\t/,$l); if($l=~m/^#/) {print "$l\tNA\n";next} if($a[0] != $a[3]) {print "$l\n";next;}  foreach $m (`cat $ARGV[1]`) {chomp($m); @b=split(/\t/,$m); next if($a[0] ne $b[0]);if($a[1]>$b[1] && $a[1]<$b[2]) {print "$l\tXXX$b[3]\n";next;} if($a[1]<$b[1] && $a[2]>$b[2]) {print "$l\tYYY$b[3]\n";next;} if($a[2]>$b[1] && $a[2]<$b[2]) {print "$l\tZZZ$b[3]\n";next;} } }' ${SV} ~/local/resources/GRCh37.e75.gene_boundaries.bed > ${NAME}.sv.bedpe 


perl -e 'foreach $l (`less $ARGV[0]`) {chomp($l); @a=split(/\t/,$l); if($l=~m/^#/) {print "$l\tNA\n";next}  foreach $m (`cat $ARGV[1]`) {chomp($m); @b=split(/\t/,$m); next if($a[0] ne $b[0]);if($a[1]>$b[1] && $a[1]<$b[2]) {print "$l\tXXX$b[3]\n";next;} if($a[1]<$b[1] && $a[2]>$b[2]) {print "$l\tYYY$b[3]\n";next;} if($a[2]>$b[1] && $a[2]<$b[2]) {print "$l\tZZZ$b[3]\n";next;} } }' ${CNV} ~/local/resources/GRCh37.e75.gene_boundaries.bed > ${NAME}.cnv.seg


done 
