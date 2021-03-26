use strict;
use warnings;


if(!(-e $ARGV[0]))
{
        die "USAGE: perl $0 INPUT_VAR\n";
}

my $rmsk="/ifs/res/leukgen/home/yellapav/resources/misc/repeat_masker.bed";
open(RMSK, "<$rmsk") or die "Can't open dbsnp file: $!\n";
open(BRASS, "<$ARGV[0]") or die "Can't open dbsnp file: $!\n";
my %repeats=();
my %samples=();
my $common=0;
while(<RMSK>) {
chomp;
next if($_ =~ m/^Hugo|^#/);
my @temp=split('\t', $_);
my $h_key=join("\t", $temp[1],$temp[2]);
unless ($repeats{$temp[0]}{$h_key}) {$repeats{$temp[0]}{$h_key}=$temp[3];}

}


while(<BRASS>) {
chomp;
my $reason1="NA";
my $reason2="NA";
if($_ =~ m/^chr1|^CHR1/) {print "$_\tRepeat_Region_BK1\tRepeat_Region_BK2\n"; next;}

my @brass=split(/\t/, $_);
   
   while (my ($key, $value) = each %{ $repeats{$brass[0]} } ) {
        chomp($key); chomp($value);
        my @temp=split(/\t/, $key);
        if(($temp[0] >= $brass[1] && $temp[0] <= $brass[2]) || ($temp[0] <= $brass[1] && $temp[1] >= $brass[2]) || ($temp[1] >= $brass[1] && $temp[1] <= $brass[2] ) || ($temp[0] >= $brass[1] && $temp[1]<= $brass[2])) {
        if($reason1 eq "NA") {$reason1=$value;} 
        else { $reason1=join(":",$reason1,$value);} 
     }
  }
   while (my ($key, $value) = each %{ $repeats{$brass[3]} } ) {
        chomp($key); chomp($value);
        my @temp=split(/\t/, $key);
        if(($temp[0] >= $brass[4] && $temp[0] <= $brass[5]) || ($temp[0] <= $brass[4] && $temp[1] >= $brass[5]) || ($temp[1] >= $brass[4] && $temp[1] <= $brass[5] ) || ($temp[0] >= $brass[4] && $temp[1]<= $brass[5])) {
        if($reason2 eq "NA") {$reason2=$value;} 
        else { $reason2=join(":",$reason2,$value);} 
     }
  }
     

print "$_\t$reason1\t$reason2\n";

} 

close(RMSK);
close(BRASS);
