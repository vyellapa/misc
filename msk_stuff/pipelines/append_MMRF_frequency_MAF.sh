MAF=$1

perl -e '%s=map{chomp;@a=split(/ /); ($a[1],$a[0])}`cat /ifs/res/leukgen/home/yellapav/resources/misc/MMRF_Canonical_GLFilter_uniq.maf|cut -f5,6,7|sort|uniq -c|sed "s/^ *//g"`; foreach $l(`cat $ARGV[0]`) {chomp($l); if($l=~m/^#/) {print "$l\n";next;} if($l=~m/^Hugo/) {print "$l\tMMRF_Frequency\n";next;} @a=split(/\t/,$l);$j=join("\t",$a[4],$a[5],$a[6]); if($s{$j}) {print "$l\t$s{$j}\n";} else {print "$l\t0\n"}}' ${MAF}
