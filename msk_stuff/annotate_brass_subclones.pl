use strict;
use warnings;


if(!(-e $ARGV[0]))
{
        die "USAGE: perl $0 INPUT_VAR\n";
}

my $rmsk="/home/yellapav/local/resources/GRCh37.e75.gene_boundaries.bed";
open(RMSK, "<$rmsk") or die "Can't open dbsnp file: $!\n";
open(BRASS, "<$ARGV[0]") or die "Can't open dbsnp file: $!\n";
my %repeats=();
my %samples=();
my $common=0;
while(<RMSK>) {
chomp;
next if($_ =~ m/^Hugo|^#/);
my @temp=split('\t', $_);
my @a = split(";",$temp[3]);
my $h_key=join("\t", $temp[1],$temp[2]);
unless ($repeats{$temp[0]}{$h_key}) {$repeats{$temp[0]}{$h_key}=$a[2];}

}


while(<BRASS>) {
chomp;
my $reason1="NA";
my $reason2="NA";
#if($_ =~ m/^chr1|^CHR1|	/) {print "$_\tAnnotation\n"; next;}

my @brass=split(/\t/, $_);
if($brass[1] =~ m/^chr|^CHR/) {print "$_\tAnnotation\n"; next;}
if(($brass[3] - $brass[2]) > 1000000) {print "$_\tNA\n"; next;}
   
   while (my ($key, $value) = each %{ $repeats{$brass[0]} } ) {
        chomp($key); chomp($value);
        my @temp=split(/\t/, $key);
        if(($temp[0] >= $brass[1] && $temp[0] <= $brass[2]) || ($temp[0] <= $brass[1] && $temp[1] >= $brass[2]) || ($temp[1] >= $brass[1] && $temp[1] <= $brass[2] ) || ($temp[0] >= $brass[1] && $temp[1]<= $brass[2])) {
        if($reason1 eq "NA") {$reason1=$value;} 
        else { $reason1=join(":",$reason1,$value);} 
     }
  }
     

print "$_\t$reason1\n";

} 

close(RMSK);
close(BRASS);
