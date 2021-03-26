perl -e 'foreach $l(`less $ARGV[0]`) {$flag=0;foreach $m(`less ~/local/msk_scripts/p219_normal_bp.txt`) {chomp($l,$m);@a=split(/\t/,$l);@b=split(/\t/,$m); if($a[0]==$b[0] and ($a[1]<$b[1]+50 and $a[1]>$b[1]-50)) {$flag=1;} if($a[3]==$b[3] and ($a[4]<($b[1]+50) and $a[4]>($b[1]-50))) {$flag=1;} } if($flag==0) {print "$l\n";} }' p162_delly_results_metrics_sep28.tsv > temp

cat temp | awk '{if(NR==1) print $0}' 
cat temp | awk '{if($1==14 && $17 > 29 && $18 > 59) print $0}'  
cat temp | awk '{if($4==14 && $14 > 29 && $15 > 59) print $0}'  
