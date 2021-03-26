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


open(VCF, "<$ARGV[0]") or die "Can't open VCF file: $!\n";
my $bam=$ARGV[1];
#open(BAMRC, "<$ARGV[1]") or die "Can't open sample file: $!\n";
my %bamrc = ();
my @temp;
my $line;
my $vcf_line;
my @REF;
my @ALT;

while(<VCF>) {

my $reason="";
chomp;
next if(($_=~/^#/) || ($_=~/^Hugo/));
my @vcf_line=split(/\t/);
my $homPosEnd=$vcf_line[1]+$homlen+1;
my $homNegEnd=$vcf_line[1]-$homlen-1;

my $hompos=join(":", $vcf_line[0],$vcf_line[1]+1);
$hompos=join("-", $hompos,$homPosEnd);

my $homneg=join(":", $vcf_line[0],$homNegEnd);
$homneg=join("-", $homneg,$vcf_line[1]-1);

my $pos=join(":", $vcf_line[0],$vcf_line[1]);
$pos=join("-", $pos,$vcf_line[1]);

my $rc=`bam-readcount -q 1 -b 10 -f /ifs/work/leukgen/ref/homo_sapiens/GRCh37d5/homo_sapiens.fasta $bam $pos`;
if($rc eq "") {warn "Seg fault for $pos ...\nMoving to next position\n";next;}

$hompos=`samtools faidx /ifs/work/leukgen/ref/homo_sapiens/GRCh37d5/homo_sapiens.fasta $hompos | grep -v "^>"`;
$homneg=`samtools faidx /ifs/work/leukgen/ref/homo_sapiens/GRCh37d5/homo_sapiens.fasta $homneg | grep -v "^>"`;
chomp($hompos); chomp($homneg);
if(maxHom($hompos)>4) {$vcf_line[6]="FAIL";$reason="HOM_POL";}
if(maxHom($homneg)>4) {$vcf_line[6]="FAIL";$reason="HOM_POL";}

@temp=split(/\t/, $rc);
my @A=split(/:/,$temp[5]);
my @C=split(/:/,$temp[6]);
my @G=split(/:/,$temp[7]);
my @T=split(/:/,$temp[8]);
my @N=split(/:/,$temp[9]);



if($vcf_line[3] eq "A") {@REF=@A;}
if($vcf_line[3] eq "C") {@REF=@C;}
if($vcf_line[3] eq "G") {@REF=@G;}
if($vcf_line[3] eq "T") {@REF=@T;}


if($vcf_line[4] eq "A") {@ALT=@A;}
if($vcf_line[4] eq "C") {@ALT=@C;}
if($vcf_line[4] eq "G") {@ALT=@G;}
if($vcf_line[4] eq "T") {@ALT=@T;}

if(abs($ALT[$avg_mapping_quality]-$REF[$avg_mapping_quality])>30) {$vcf_line[6]="FAIL";$reason=join(":",$reason,"MAP_QUAL_DIFF");}
if(abs($ALT[$avg_mapping_quality])<29) {$vcf_line[6]="FAIL";$reason=join(":",$reason,"MAP_QUAL_ALT");}
if(abs($REF[$avg_mapping_quality])<29) {$vcf_line[6]="FAIL";$reason=join(":",$reason,"MAP_QUAL_REF");}
if(abs($ALT[$avg_sum_mismatch_qualities]-$REF[$avg_sum_mismatch_qualities])>50) {$vcf_line[6]="FAIL";$reason=join(":",$reason,"MMQS_DIFF");}
if(abs($REF[$avg_sum_mismatch_qualities])>60) {$vcf_line[6]="FAIL";$reason=join(":",$reason,"MMQS_REF");}
if(abs($ALT[$avg_sum_mismatch_qualities])>60) {$vcf_line[6]="FAIL";$reason=join(":",$reason,"MMQS_ALT");}
if(abs($REF[$count])<3) {$vcf_line[6]="FAIL";$reason=join(":",$reason,"ALT_COUNT");}
#if(abs($ALT[$avg_mapping_quality]-$REF[$avg_mapping_quality])>29) {$reason=join(":",$reason,"MAP_QUAL_DIFF");}

chomp($hompos); chomp($homneg); chomp($rc);
my $line=join("\t",@vcf_line);
print "$line\t$reason\t$rc\t$homneg\t$hompos\n";

}

close(VCF);




sub maxHom {

my $seq=shift;
$seq=~ tr/acgt/ACGT/;
my $ap=$seq=~tr/A//;
my $ac=$seq=~tr/C//;
my $ag=$seq=~tr/G//;
my $at=$seq=~tr/T//;
my $maxn=max($ap,$ac,$ag,$at);
return $maxn;

}

