
UD=$1
TSV=$2



perl -e '%s=map{chomp;@a=split(/\t/);($a[0],1)}`less $ARGV[0]`; foreach $l (`less $ARGV[1]`) {chomp($l); @a=split(/\t/,$l); if($s{$a[0]}) {print "$l\n";}}' ${UD} ${TSV} > t

cat t| awk -F'\t' '{if($22!~/^synon/ && $22!~/intron/ && $22!~/stream/ && $22!~/prime/ && $22!="NA") print $0}'
