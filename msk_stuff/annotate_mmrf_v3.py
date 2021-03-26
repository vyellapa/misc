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
def annotate_tsv_freq(in_tsv_gz,annotation_tsv):
    """Reading and creating vcf files."""
    warnings.warn("Reading TSV file ...")
    #print(counts_tsv.head())


    mmrf = pd.read_csv('/ifs/res/leukgen/home/yellapav/MMRF/MMRF_CoMMpass_IA9_All_Canonical_Variants.txt', sep="\t")
    mmrf=mmrf.iloc[:,[0,1,2,4,5,19,23]]
    mmrf=mmrf.drop_duplicates()

    mmrfM=mmrf.groupby(['CHROM','POS'])['GEN[1].AR'].median()
    mmrfC=mmrf.groupby(['CHROM','POS'])['GEN[1].AR'].count()
    mmrfQ25=mmrf.groupby(['CHROM','POS'])['GEN[1].AR'].quantile(q=0.25)
    mmrfQ75=mmrf.groupby(['CHROM','POS'])['GEN[1].AR'].quantile(q=0.75)
    

    #anno_tsv = pd.read_csv(annotation_tsv, comment='#',sep="\t")
    anno_tsv = pd.read_csv(annotation_tsv, comment='#',sep="\t", low_memory=False)
    #anno_tsv[anno_tsv['FILTER'] == "PASS"]
    counts_tsv=anno_tsv.groupby(["CHR","START","REF","ALT"]).size().reset_index(name="count")
    counts_tsv=counts_tsv[["CHR", "START","count"]].set_index(['CHR','START'])
    counts_median=anno_tsv.groupby(['CHR','START'])['TARGET_VAF'].median()









    inFile = gzip.open(in_tsv_gz,'r')
    
    warnings.warn("Annotating ...")
    for record in inFile:
        record=record.decode("utf-8")
        record=record.rstrip()
        recArr=record.split("\t")
        
        cl = [] 
        freq = [] 
        medVAF = [] 
        Q25 = [] 
        Q75 = [] 
        positions = [] 
        normal = "0" 
        normalVAF = "0" 
        if record.startswith("#"):
          continue

        if recArr[0] == "ID_VARIANT":
           cl = "MMRF_Class"
           freq = "MMRF_Frequency"
           medVAF = "MMRF_VAF"
           Q25 = "MMRF_Q25"
           Q75 = "MMRF_Q75"
           positions = "MMRF_Positions"
           normal = "Normals_Frequency"
           normalVAF = "Normals_median_VAF"
           record = [ record, cl, freq, medVAF, Q25, Q75, positions, normal, normalVAF ]
           record = ("\t".join(record))
           print(record)
           continue

        try:
            chrom = str(recArr[3])
            pos = int(recArr[4])
            start = int(recArr[4]) - 9
            end = int(recArr[4]) + 9
            if (chrom, pos) in mmrfC.index:
                cl = "genomic_exact"
                freq = str(mmrfC.loc[(chrom,pos)]) 
                medVAF = str(mmrfM.loc[(chrom,pos)]) 
                Q25 = str(mmrfQ25.loc[(chrom,pos)]) 
                Q75 = str(mmrfQ75.loc[(chrom,pos)]) 
                positions = str(pos)
                record = [ record, cl, freq, medVAF, Q25, Q75, positions ]
                record = ("\t".join(record))
                #print(record)
                flag = 1
            if flag == 0: 
                mmrfCsub=mmrfC.loc[chrom]
                if not mmrfCsub[(mmrfCsub.index >= start) & (mmrfCsub.index <= end)].empty:
                    for i in mmrfCsub[(mmrfCsub.index >= start) & (mmrfCsub.index <= end)].index.values:
                        #print(i, mmrfC.loc[(chrom,i)], start, end)
                        cl = "genomic_close"
                        freq.append(str(mmrfC.loc[(chrom,i)]))
                        medVAF.append(str(mmrfM.loc[(chrom,i)]))
                        Q25.append(str(mmrfQ25.loc[(chrom,i)]))
                        Q75.append(str(mmrfQ75.loc[(chrom,i)]))
                        positions.append(str(i))
                    freq = (":".join(freq))
                    medVAF = (":".join(medVAF))
                    Q25 = (":".join(Q25))
                    Q75 = (":".join(Q75))
                    positions = (":".join(positions))
                    record = [ record, cl, freq, medVAF, Q25, Q75, positions ]
                    record = ("\t".join(record))
                    #print(record)
                    #continue
                else:
                    cl = "NA"
                    freq = "NA"
                    medVAF = "NA"
                    Q25 = "NA"
                    Q75 = "NA"
                    positions = "NA"
                    record = [ record, cl, freq, medVAF, Q25, Q75, positions ]
                    record = ("\t".join(record))
                    #print(record)
                    #continue


        except:
            cl = "NA"
            freq = "NA"
            medVAF = "NA"
            Q25 = "NA"
            Q75 = "NA"
            positions = "NA"
            record = [ record, cl, freq, medVAF, Q25, Q75, positions ]
            record = ("\t".join(record))
            #print(record)


        normal = "0"
        normalVAF = "0"
        try:
            chrom=str(recArr[3])
            pos=int(recArr[4])
            #print(chrom)
            #print(pos)
            normal = counts_tsv.loc[(chrom,pos),"count"]
            normal = normal.ix[0]
            normal = str(normal)

            normalVAF = str(counts_median.loc[(chrom,pos)])

            record = [ record, normal, normalVAF ]
            record = ("\t".join(record))
            #print(type(freq))
            print(record)
            #print(freq.iloc[0:1,0:1])

        except:
            normal = "0"
            normalVAF = "0"
            record = [ record, str(normal), str(normalVAF) ]
            record = ("\t".join(record))
            print(record)





@click.group()
def commands():
    """Usage for TSV annotator."""
    pass


commands.add_command(annotate_tsv_freq)

if __name__ == '__main__':
    commands()
