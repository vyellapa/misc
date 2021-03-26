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
BOLLI = '/home/yellapav/local/resources/bolli_mutations.txt'
### Add RG Tags
@click.command()
@click.option('--in_tsv_gz', required=True, help='Input TSV File')
@click.option('--annotation_tsv', required=True, help='Annotation TSV File')
def annotate_tsv_freq(in_tsv_gz,annotation_tsv):
    """Reading and creating vcf files."""
    sys.stderr.write("Reading TSV file ...\n")
    nicollo  = pd.read_csv(BOLLI, sep="\t")
    nicollo = nicollo.iloc[:,[1,2,4,5,23]]
    nicollo_counts = nicollo.groupby(['CHR','START'])['MT'].count()
    nol_var = nicollo.drop(['WT','MT'], axis = 1) 
    nol_var = nol_var.set_index(['CHR', 'START'])

    #nicollo_counts = nicollo.groupby(["CHR","START","WT","MT"]).size().reset_index(name="count")
    #nicollo_counts = nicollo_counts[["CHR", "START","count"]].set_index(['CHR','START'])

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
    
    sys.stderr.write("Annotating ...\n")
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
        bolli_cl = [] 
        bolli_freq = [] 
        bolli_positions = [] 
        bolli_anno = [] 
        flag = 0
        bolli_flag = 0
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
           bolli_cl = "Bolli_Class"
           bolli_freq = "Bolli_Frequency"
           bolli_positions = "Bolli_Positions"
           bolli_anno = "Bolli_Annotation"
           record = [ record, cl, freq, medVAF, Q25, Q75, positions, bolli_cl, bolli_freq, bolli_anno, bolli_positions, normal, normalVAF ]
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
                flag = 1
            if flag == 0:
                mmrfCsub=mmrfC.loc[chrom]
                if not mmrfCsub[(mmrfCsub.index >= start) & (mmrfCsub.index <= end)].empty:
                    for i in mmrfCsub[(mmrfCsub.index >= start) & (mmrfCsub.index <= end)].index.values:
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
                else:
                    cl = "NA"
                    freq = "NA"
                    medVAF = "NA"
                    Q25 = "NA"
                    Q75 = "NA"
                    positions = "NA"
                    record = [ record, cl, freq, medVAF, Q25, Q75, positions ]
                    record = ("\t".join(record))


        except:
            cl = "NA"
            freq = "NA"
            medVAF = "NA"
            Q25 = "NA"
            Q75 = "NA"
            positions = "NA"
            record = [ record, cl, freq, medVAF, Q25, Q75, positions ]
            record = ("\t".join(record))



        try:
            chrom = str(recArr[3])
            pos = int(recArr[4])
            start = int(recArr[4]) - 9
            end = int(recArr[4]) + 9
            if (chrom, pos) in  nicollo_counts.index:
                bolli_cl = "genomic_exact"
                bolli_freq = str(nicollo_counts.loc[(chrom,pos)]) 
                bolli_positions = str(pos)
                bolli_anno = str(nol_var.loc[chrom, pos]['Variant_class'].values[0])
                record = [ record, bolli_cl, bolli_freq, bolli_anno, bolli_positions ]
                record = ("\t".join(record))
                bolli_flag = 1


            if bolli_flag == 0: 
                nicollo_counts_sub=nicollo_counts.loc[chrom]
                if not nicollo_counts_sub[(nicollo_counts_sub.index >= start) & (nicollo_counts_sub.index <= end)].empty:
                    for i in nicollo_counts_sub[(nicollo_counts_sub.index >= start) & (nicollo_counts_sub.index <= end)].index.values:
                #if not nicollo_counts_sub.ix[start:end].empty:
                #    for i in nicollo_counts_sub.ix[start:end].index.values:
                        #print("XXXXXXX",i, nicollo_counts_sub.loc[(chrom,i)], start, end)
                        bolli_cl = "genomic_close"
                        bolli_freq.append(str(nicollo_counts.loc[(chrom,i)]))
                        bolli_anno.append(str(nol_var.loc[(chrom,i)]['Variant_class'].values[0]))
                        bolli_positions.append(str(i))
                    bolli_freq = (":".join(bolli_freq))
                    bolli_positions = (":".join(bolli_positions))
                    bolli_anno = (":".join(bolli_anno))
                    record = [ record, bolli_cl, bolli_freq, bolli_anno, bolli_positions ]
                    record = ("\t".join(record))
                else:
                    bolli_cl = "NA"
                    bolli_freq = "NA"
                    bolli_positions = "NA"
                    bolli_anno = "NA"
                    record = [ record, bolli_cl, bolli_freq, bolli_anno, bolli_positions ]
                    record = ("\t".join(record))


        except:
            bolli_cl = "NA"
            bolli_freq = "NA"
            bolli_anno = "NA"
            bolli_positions = "NA"
            record = [ record, bolli_cl, bolli_freq, bolli_anno, bolli_positions ]
            record = ("\t".join(record))


        normal = "0"
        normalVAF = "0"
        try:
            chrom=str(recArr[3])
            pos=int(recArr[4])
            normal = counts_tsv.loc[(chrom,pos),"count"]
            normal = normal.ix[0]
            normal = str(normal)

            normalVAF = str(counts_median.loc[(chrom,pos)])

            record = [ record, normal, normalVAF ]
            record = ("\t".join(record))
            print(record)

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
