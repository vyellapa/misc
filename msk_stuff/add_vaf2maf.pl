#! /gsc/bin/perl
################################################################
## exclude.pl - Takes in two mafs and excludes                ##
## missense mutations that are below a certain threshold      ##
##      AUTHOR:vyellapa                                       ##
##      Last updated:Fri Oct 17 17:35:33 CDT 2014             ##
################################################################

use strict;
use warnings;
use FileHandle;
use File::Basename;
use List::Util qw[min max];

if(!(-e $ARGV[1]))
{
        die "USAGE: perl $0 VCF BAM\n";
}

my $homlen=5;
my $base=0;
my $count=1;
my $avg_mapping_quality=2;
my $avg_basequality=3;
my $avg_se_mapping_quality=4;
my $num_plus_strand=5;
my $num_minus_strand=6;
my $avg_pos_as_fraction=7;
my $avg_num_mismatches_as_fraction=8;
my $avg_sum_mismatch_qualities=9;
my $num_q2_containing_reads=10;
my $avg_distance_to_q2_start_in_q2_reads=11;
my $avg_clipped_length=12;
my $avg_distance_to_effective_3p_end=13;


open(MAF, "<$ARGV[0]") or die "Can't open VCF file: $!\n";
my $bam=$ARGV[1];
my $name=$ARGV[2];
#open(BAMRC, "<$ARGV[1]") or die "Can't open sample file: $!\n";
my %bamrc = ();
my @temp;
my $line;
my $maf_line;
my @REF;
my @ALT;

while(<MAF>) {

my $vaf="";
chomp;
my @maf_line=split(/\t/);
my $size=scalar(@maf_line);
if(($_=~/^#/) || ($_=~/^Hugo/)) {print "$_\t$name\n";next;}

my $pos=join(":", $maf_line[4],$maf_line[5]);
$pos=join("-", $pos,$maf_line[6]);

my $rc=`bam-readcount -q 30 -b 15 -f /ifs/work/leukgen/ref/homo_sapiens/GRCh37d5/genome/gr37.fasta $bam $pos`;
if($rc eq "") {warn "Seg fault for $pos ...\nMoving to next position\n";next;}


@temp=split(/\t/, $rc);
my @A=split(/:/,$temp[5]);
my @C=split(/:/,$temp[6]);
my @G=split(/:/,$temp[7]);
my @T=split(/:/,$temp[8]);
my @N=split(/:/,$temp[9]);

if($maf_line[10] eq "A") {@REF=@A;}
if($maf_line[10] eq "C") {@REF=@C;}
if($maf_line[10] eq "G") {@REF=@G;}
if($maf_line[10] eq "T") {@REF=@T;}

if($maf_line[12] eq "A") {@ALT=@A;}
if($maf_line[12] eq "C") {@ALT=@C;}
if($maf_line[12] eq "G") {@ALT=@G;}
if($maf_line[12] eq "T") {@ALT=@T;}

$maf_line[$size]=0;
if($ALT[$count]>0 || $REF[$count]>0) {
$maf_line[$size]=($ALT[$count]/($ALT[$count]+$REF[$count]));}
#if(abs($ALT[$avg_mapping_quality]-$REF[$avg_mapping_quality])>29) {$reason=join(":",$reason,"MAP_QUAL_DIFF");}

my $line=join("\t",@maf_line);
print "$line\n";

}

close(MAF);



