tsv=$1
vcf=$2

perl -e '%s=map{chomp;@a=split(/\t/);($a[2],1)}`less $ARGV[0]|grep -v "^#"`;$col=2;if($ARGV[1]=~m/\.maf/) {$col=110;} foreach $l (`less $ARGV[1]`) {chomp($l); @a=split(/\t/,$l);if($l=~m/^ID_VARIANT/) {print "$l\n"; next;} if($s{$a[0]}) {print "$l\n";} }' ${vcf} ${tsv}
