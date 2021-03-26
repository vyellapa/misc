#perl brass_repeatMask.pl /ifs/res/leukgen/local/opt/leukdc/data/analyses/02/56/60256/output/E-H-131564-T1-1-D1-1.cns
use strict;
use warnings;


if(!(-e $ARGV[0]))
{
        die "USAGE: perl $0 INPUT_VAR\n";
}

#my $rmsk="/ifs/res/leukgen/home/yellapav/resources/misc/repeat_masker.bed";
my $rmsk="/home/yellapav/local/resources/GRCh37.e75.gene_boundaries.bed";
open(RMSK, "<$rmsk") or die "Can't open dbsnp file: $!\n";
open(BRASS, "<$ARGV[0]") or die "Can't open dbsnp file: $!\n";
my %repeats=();
my %samples=();
my %final=();
my $common=0;
while(<RMSK>) {
chomp;
next if($_ =~ m/^Hugo|^#/);
my @temp=split('\t', $_);
my @gene=split(';', $temp[3]);
my $h_key=join("\t", $temp[1],$temp[2]);
unless ($repeats{$temp[0]}{$h_key}) {$repeats{$temp[0]}{$h_key}=$gene[2];}

}


while(<BRASS>) {
chomp;
my $reason1="NA";
my $reason2="NA";
if($_ =~ m/^chromos|^CHR1/) {print "Gene\tCN\n"; next;}

my @brass=split(/\t/, $_);
   
   while (my ($key, $value) = each %{ $repeats{$brass[0]} } ) {
        chomp($key); chomp($value);
        my @temp=split(/\t/, $key);
        if(($temp[0] >= $brass[1] && $temp[0] <= $brass[2]) || ($temp[0] <= $brass[1] && $temp[1] >= $brass[2]) || ($temp[1] >= $brass[1] && $temp[1] <= $brass[2] ) || ($temp[0] >= $brass[1] && $temp[1]<= $brass[2])) {
              my $cn=2*(2**$brass[4]);
              unless($final{$value}) {$final{$value}=$cn;}
              if($final{$value} < $cn && $final{$value}) {$final{$value} = $cn;}
              #print "$value\t$cn\t$brass[4]\n";
     }
  }

}

foreach my $gene (`cat $rmsk |cut -d';' -f3|sort -u`) 
{
chomp($gene);
if($final{$gene}) {print "$gene\t$final{$gene}\n";}
else {print "$gene\t2\n";}



}




close(RMSK);
close(BRASS);
