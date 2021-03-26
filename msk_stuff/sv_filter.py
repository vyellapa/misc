import pysam
import re
import click

@click.command()
@click.option(
    "--tumor_bam",
    show_default=True,
    type=click.Path(exists=True),
    help="Path to tumour bam file.",
    default="/ifs/res/leukgen/local/opt/leukdc/data/workflows/31/04/33104/data/bam/I-H-130720-T1-1-D1-1.bam"
    # required=True,
    )

@click.option(
    "--bcf",
    show_default=True,
    type=click.Path(exists=True),
    help="Path to Delly translocation bcf file.",
    default="/ifs/res/leukgen/projects/220/RESULTS/SV/delly_220/I-H-130720-T1-1-D1-1_TRA.bcf"
    # required=True,
    )


@click.option(
    "--panel_of_normals",
    show_default=True,
    type=click.Path(exists=True),
    help="Path to Delly panel of normals FP break points tsv file.(CHR POS)",
    default="/ifs/res/leukgen/projects/220/RESULTS/SV/delly_220/I-H-130720-T1-1-D1-1_TRA.bcf"
    # required=True,
    )

@click.option(
    "--gtf",
    show_default=True,
    type=click.Path(exists=True),
    help="Path to the gtf for annotation.\n This should be tabixed.",
    default="/home/yellapav/local/resources/Homo_sapiens.GRCh37.75.sorted.gtf.gz"
    # required=True,
    )

@click.option(
    "--sample",
    show_default=True,
    help="Sample name",
    default="I-H-130720-T1-1-D1-1"
    # required=True,
    )

@click.option(
    "--spanning_pairs",
    show_default=True,
    type=int,
    help="Number of spanning pairs below which SVs will be filtered.",
    default=5
    )

@click.option(
    "--junction_reads",
    show_default=True,
    type=int,
    help="Number of junction reads below which SVs will be filtered.",
    default=1
    )

@click.option(
    "--spanning_pairs_only",
    show_default=True,
    type=int,
    help="Number of spanning pairs to rescue SVs with regardless of thier junction reads.",
    default=30
    )

@click.option(
    "--fragment_size",
    show_default=True,
    type=int,
    help="Fragment Size",
    default=500
    )


def filter_delly(tumor_bam,bcf,panel_of_normals,gtf,sample,spanning_pairs,junction_reads,spanning_pairs_only,fragment_size):
    
    tabixfile = pysam.TabixFile(gtf)
    fragment_size = int(fragment_size) * 2


    print("Chrom1","Pos1","End1","Chrom2","Pos2","End2","SV_Type","Strand1","Strand2","Sample","Tumor_Juntion_Pairs", "Tumor_Spanning_Pairs","Normal_Juntion_Pairs", "Normal_Spanning_Pairs","Support", "MAPQ", "MATCH", "Annotation_bp1", "Annotation_bp2",sep="\t")


    bcf_in = pysam.VariantFile(bcf)
    for rec in bcf_in.fetch():
        flag1 = 0
        flag2 = 0
        flag3 = 0
        flag4 = 0
        support = 0
        match = 0
        mapq = 0
        tumor = 0
        normal = 1
        strand1 = "NA"
        strand2 = "NA"
        sv_type = "translocation"
        annotation_bp1 = []
        annotation_bp2 = []
        annotation_bp1_str = "NA"
        annotation_bp2_str = "NA"
        i = 0
        #filter if not PASS
        filter = ''.join(list(rec.filter))
        if filter != "PASS":
            continue

    # filter if either contig is non 14
        if rec.contig != "14" and rec.info['CHR2'] != "14":
            continue

    # Skip if it does not pass the spanning reads and junctions 
        if rec.samples[tumor]["DV"] <= spanning_pairs and rec.samples[0]["RV"] <= junction_reads:
            flag1 = 1

        if rec.samples[tumor]["DV"] <= spanning_pairs_only:
            flag2 = 1

        if flag1 == 1 and flag2 == 1:
            continue

        if (rec.samples[tumor]["GT"]) == (rec.samples[normal]["GT"]):
            #print(rec)
            continue

        with open(panel_of_normals, 'r', encoding="ISO-8859-1") as pon:
            for line in pon:
                line = line.rstrip()
                line = re.split(r'\t+', line)
                if rec.info['CHR2'] == line[0] and rec.stop > (int(line[1])-50) and rec.stop < (int(line[1])-50):
                    flag3 = 1

                if rec.contig == line[0] and rec.pos > (int(line[1])-50) and rec.pos < (int(line[1])-50):
                    flag4 = 1

        if flag3 == 1 or flag4 == 1:
            continue
        
        if rec.info['CT'] == "3to5":
            strand1 = "+"
            strand2 = "+"
        
        if rec.info['CT'] == "5to3":
            strand1 = "-"
            strand2 = "-"
        
        if rec.info['CT'] == "5to5":
            strand1 = "-"
            strand2 = "+"
        
        if rec.info['CT'] == "3to3":
            strand1 = "+"
            strand2 = "-"

        samfile = pysam.AlignmentFile(tumor_bam, "rb")

        if rec.contig == "14":
            igh_partner_chrom = rec.info['CHR2']
            igh_partner_pos = rec.stop
            mate_pos = rec.pos
        else:
            igh_partner_chrom = rec.contig
            igh_partner_pos = rec.pos
            mate_pos = rec.stop

        for read in samfile.fetch(igh_partner_chrom, igh_partner_pos - fragment_size, igh_partner_pos + fragment_size):
            mate = samfile.mate(read)
            if mate.reference_name != "14":
                continue

            if int(mate.reference_start) > int(mate_pos)-fragment_size and int(mate.reference_start) < int(mate_pos)+fragment_size:
                i = i + 1
                #print("YES",mate.reference_start, mate_pos)
            else:
                continue
	#if
            support += 1
            mapq += read.mapping_quality
            match_cigar = re.findall('(\d+)M', read.cigarstring)
            match += sum(map(int, match_cigar))

        #print(read.cigarstring, support)
        #print(read.mapping_quality, mate.reference_start)
        #print("SEE",support, mapq, match, match)

        samfile.close()
        mapq = mapq/support
        match = match/support
    
        for row in tabixfile.fetch(rec.contig, int(rec.pos)-1, rec.pos, parser=pysam.asGTF()):
            annotation_bp1.append(row.gene_name + ";" + row.gene_id + ";" + row.strand)
            annotation_bp1_str = "|".join((list(set(annotation_bp1))))
        
        for row in tabixfile.fetch(rec.info['CHR2'], int(rec.stop)-1, rec.stop, parser=pysam.asGTF()):
            annotation_bp2.append(row.gene_name + ";" + row.gene_id + ";" + row.strand)
            annotation_bp2_str = "|".join((list(set(annotation_bp2))))
        
    #print(sample,rec.contig,rec.pos,rec.info['CHR2'],rec.stop,rec.samples[tumor]["DV"], rec.samples[tumor]["RV"], rec.samples[normal]["DV"], rec.samples[normal]["RV"], support, mapq, match, annotation_bp1_str, annotation_bp2_str,sep="\t")
        print(rec.contig,rec.pos,int(rec.pos)+100,rec.info['CHR2'],rec.stop,int(rec.stop)+100,sv_type,strand1,strand2,sample,rec.samples[tumor]["DV"], rec.samples[tumor]["RV"], rec.samples[normal]["DV"], rec.samples[normal]["RV"], support, mapq, match, annotation_bp1_str, annotation_bp2_str,sep="\t")
    
    
    
if __name__ == '__main__':
    filter_delly()
