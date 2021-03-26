"""Usage python annotate_delly.py annotate_sv_vcf --in_vcf ~/example.vcf --out_vcf hey.vcf"""

import sys
import os
import subprocess
import click
import vcf

GTF = '/home/yellapav/local/resources/Homo_sapiens.GRCh37.75.gtf'
BED = '/home/yellapav/local/resources/GRCh37.e75.gene_boundaries.bed'
TEMP = '/ifs/work/leukgen/home/yellapav/star_align/tmp'


f = open('/home/yellapav/local/resources/GRCh37.e75.gene_boundaries.bed','r')
print("Reading bed file ...")
gene = {}
for line in f:
    line=line.rstrip()
    arr=line.split("\t")
    gene[arr[0]]={}

f.seek(0)

for line in f:
    line=line.rstrip()
    arr=line.split("\t")
    gene[arr[0]][arr[3]]=line



### Add RG Tags
@click.command()
@click.option('--in_vcf', required=True, help='Input Vcf File')
@click.option('--out_vcf', required=True, help='Output Vcf File')
def annotate_sv_vcf(in_vcf,out_vcf):
    """Reading and creating vcf files."""
    print("Reading Vcf file ...")
    vcf_reader=vcf.Reader(open(in_vcf,'r'))
    vcf_writer=vcf.Writer(open(out_vcf,'w'), vcf_reader)
    print("Annotating ...")
    for record in vcf_reader:
        #record=record.rstrip()

        for key,value in gene[record.CHROM].items():
            arr=value.split("\t")
            if int(arr[1])<=record.POS and int(arr[2])>=record.POS:
                new = [ record.INFO['SVMETHOD'],"BK1", arr[3] ]
                new=(",".join(new))
                new=new.replace(";",",")
                record.INFO['SVMETHOD']=new

        for key,value in gene[record.INFO['CHR2']].items():
            arr=value.split("\t")
            if int(arr[1])<=record.INFO['END'] and int(arr[2])>=record.INFO['END']:
                new = [ record.INFO['SVMETHOD'],"BK2", arr[3] ]
                new=(",".join(new))
                new=new.replace(";",",")
                record.INFO['SVMETHOD']=new


        #print(record.INFO['SVMETHOD'])
        vcf_writer.write_record(record)
    vcf_writer.close()




@click.group()
def commands():
    """USage for SV annotator."""
    pass


commands.add_command(annotate_sv_vcf)

if __name__ == '__main__':
    commands()
