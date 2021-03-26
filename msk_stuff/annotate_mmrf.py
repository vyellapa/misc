"""Usage python annotate_tsv.py annotate_with_mmrf --in_vcf ~/example.vcf --out_vcf hey.vcf"""

import sys
import os
import subprocess
import click
import vcf
import pandas as pd
import gzip 


GTF = '/home/yellapav/local/resources/Homo_sapiens.GRCh37.75.gtf'
BED = '/home/yellapav/local/resources/GRCh37.e75.gene_boundaries.bed'
TEMP = '/ifs/work/leukgen/home/yellapav/star_align/tmp'
MMRF_TSV = '/home/yellapav/local/resources/MMRF_annotateCavemans.tsv'
MMRF_TSV = '/home/yellapav/local/resources/MMRF_annotateCavemans_freq.tsv'


t = pd.read_csv('/home/yellapav/local/resources/MMRF_annotateCavemans_freq.tsv', sep = "\t")
f = open(MMRF_TSV,'r')
print("Reading MMRF mutations ...")
freq = {}
for line in f:
    line=line.rstrip()
    arr=line.split(" ")
    freq[arr[1]]={}

f.seek(0)

for line in f:
    line=line.rstrip()
    arr=line.split(" ")
    freq[arr[1]][arr[2]]=arr[0]



### Add RG Tags
@click.command()
@click.option('--in_tsv_gz', required=True, help='Input tsv file')
#@click.option('--out_vcf', required=True, help='Output Vcf File')
def annotate_with_mmrf(in_tsv_gz):
    """Reading and annotating tsv file."""
    print("Reading input file ...")
    inFile = gzip.open(in_tsv_gz,'r')
    
    print("Annotating ...")
    for record in inFile:
        record=record.decode("utf-8")
        record=record.rstrip()
        recArr=record.split("\t")
           
        mmrf = "None"
        if record.startswith("#"):
          continue
        if recArr[0] == "ID_VARIANT":
           mmrf="MMRF_Frequency"
        if recArr[3] in freq:
         for key,value in freq[recArr[3]].items():
            arr=value.split("\t")
            if int(key) == int(recArr[4]):
                mmrf = [ "genomic_exact", freq[recArr[3]][recArr[4]] ]
                mmrf = (":".join(mmrf))
                break

            ## Genomic Close
            start = int(key) - 9
            end = int(key) + 9
            
            if int(start) <= int(recArr[4]) and int(end) >= int(recArr[4]):
                if mmrf == "None":
                    mmrf = [ "genomic_close", freq[recArr[3]][key] ]
                    mmrf = (":".join(mmrf))
                else:
                    append = [ mmrf , ";genomic_close", freq[recArr[3]][key] ]
                    mmrf = (":".join(append))
                    mmrf = mmrf.replace(":;",";")
    
        record = [ record, mmrf ]
        record = ("\t".join(record))
        print(record)



@click.group()
def commands():
    """Usage for mmrf annotator."""
    pass


commands.add_command(annotate_with_mmrf)


if __name__ == '__main__':
    commands()
