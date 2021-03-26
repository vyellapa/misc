cat alll |awk -F'\t' '{OFS="\t";if($1==14 || $3==14 || NR==1) {print $0}}' | awk -F'\t' '{flag=0; if(NR==1) {print $0;flag=1} if($1==14 && ($28<29 || $29<60 || $30<4)) {flag=1;} if($3==14 && ($25<29 || $26<60 || $27<4)) {flag=1;} if(flag==0) {print $0}}' > temp

perl -e 'foreach $l(`less $ARGV[0]`) {$flag=0;foreach $m(`less ~/local/msk_scripts/brass_FP_coords.txt`) {chomp($l,$m);@a=split(/\t/,$l);@b=split(/\t/,$m); if($a[0]==$b[0] and ($a[1]<$b[1]+100 and $a[1]>$b[1]-100)) {$flag=1;} if($a[3]==$b[0] and ($a[4]<($b[1]+100) and $a[4]>($b[1]-100))) {$flag=1;} } if($flag==0) {print "$l\n";} }' temp




less p222_brassR3_50bp_results_feb8_metrics_filter.tsv | awk -F'\t' '{print "t("$1":"$3")\t"$9}' > temp


perl -e '%s={}; foreach $l(`cat temp`) {chomp($l);@a=split(/\t/,$l);$c=substr($a[1],0,15);$s{$a[0]}{$c}=1; } @a=qw(t(11:14) t(14:16) t(14:20) t(4:14) t(6:14) t(8:14)); foreach $m (`less $ARGV[0]`) {chomp($m);@z=split(/\t/,$m); if($m=~/^Al/) {$j=join("\t",@a); print "$m\tBRASS\t$j\n";next;} print "$m";foreach $l (@a) {if($s{$l}{$z[0]}) {print "\t1";} else {print "\t0";} } print "\n";}' ../../../delly/malin_results_seq.txt|cut -f109-|less

