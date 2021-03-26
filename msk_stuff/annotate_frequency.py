"""Usage python annotate_delly.py annotate_sv_vcf --in_vcf ~/example.vcf --out_vcf hey.vcf"""

import sys
import os
import subprocess
import click
import vcf
import pandas as pd
import gzip
import warnings

GTF = '/home/yellapav/local/resources/Homo_sapiens.GRCh37.75.gtf'
BED = '/home/yellapav/local/resources/GRCh37.e75.gene_boundaries.bed'
TEMP = '/ifs/work/leukgen/home/yellapav/star_align/tmp'



### Add RG Tags
@click.command()
@click.option('--in_tsv_gz', required=True, help='Input TSV File')
@click.option('--annotation_tsv', required=True, help='Annotation TSV File')
def annotate_tsv_freq(annotation_tsv,in_tsv_gz):
    """Reading and creating vcf files."""
    warnings.warn("Reading Vcf file ...")
    anno_tsv = pd.read_csv(annotation_tsv, comment='#',sep="\t")
    #anno_tsv[anno_tsv['FILTER'] == "PASS"]
    counts_tsv=anno_tsv.groupby(["CHR","START","REF","ALT"]).size().reset_index(name="count")
    counts_tsv=counts_tsv[["CHR", "START","count"]].set_index(['CHR','START'])
    #print(counts_tsv.head())
    inFile = gzip.open(in_tsv_gz,'r')
    
    warnings.warn("Annotating ...")
    for record in inFile:
        record=record.decode("utf-8")
        record=record.rstrip()
        recArr=record.split("\t")
           
        freq = "None"
        if record.startswith("#"):
          continue

        if recArr[0] == "ID_VARIANT":
           freq = "Normals_Frequency"
           record = [ record, freq ]
           record = ("\t".join(record))
           print(record)
           continue

        try:
            chrom=str(recArr[3])
            pos=int(recArr[4])
            #print(chrom)
            #print(pos)
            freq = counts_tsv.loc[(chrom,pos),"count"]
            freq = freq.ix[0]
            freq = str(freq)
            record = [ record, freq ]
            record = ("\t".join(record))
            #print(type(freq))
            print(record)
            #print(freq.iloc[0:1,0:1])

        except:
            freq = "0"
            record = [ record, str(freq) ]
            record = ("\t".join(record))
            print(record)




@click.group()
def commands():
    """Usage for TSV annotator."""
    pass


commands.add_command(annotate_tsv_freq)

if __name__ == '__main__':
    commands()
