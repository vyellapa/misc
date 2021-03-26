ids=$1
vcf=$2

perl -e '%s=map{chomp;@a=split(/\t/);($a[0],1)}`less $ARGV[0]`;$col=2;if($ARGV[1]=~m/\.maf/) {$col=110;} foreach $l (`less $ARGV[1]`) {chomp($l); @a=split(/\t/,$l);if($l=~m/^Hugo|^#/) {print "$l\n"; next;} if($s{$a[$col]}) {print "$l\n";} }' ${ids} ${vcf}
