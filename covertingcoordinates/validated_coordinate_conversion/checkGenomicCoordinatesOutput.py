#!/bin/python
import pybedtools
from pybedtools import BedTool

# First load in hg38 fasta file
fasta = pybedtools.example_filename("/data/data1/GenomeReferences/Homo_sapiens.GRCh38.dna.toplevel.fa")


import pandas as pd

# Get the final DF file to verify
dataToVerify = pd.read_csv("../tmp/openASO_gr.tsv",sep="\t",header=0)

# Loop through every row
for index, row in x.iterrows():
    chrm = str(row["chrom"])
    #Make it like a bed file
    start = str(row["chromStart"]-1)
    end = str(row["chromEnd"])
    strand = row["strand"]
    
    query = " ".join([chrm,start,end])
    a = pybedtools.BedTool(query, from_string=True)
    a = a.sequence(fi=fasta)
    seqToCheck = open(a.seqfn).read().strip().split("\n")[1]
    
    aso_seq = row["ASOseq"]
    aso_seq_rc = row["aso_rc_seq"]
    # seqToCheck is always going to be on the plus strand
    # If gene is on negative strand then the aso is going to on the plus strand
    # If gene is on positive strand then ASO seq will be on the negative strand, so we have to check the reverse complement of that [which is on the positive strand]
    if (strand=="-" and seqToCheck != aso_seq) or (strand=="+"and seqToCheck != aso_seq_rc):
        print(query)
        print(row["uniq_id"]+"\t"+seqToCheck+"\t"+aso_seq_rc+"\t"+aso_seq)
    