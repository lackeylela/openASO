# checks for genes that don't have a single aso match

# import necessary modules
import pandas as pd
from Bio.Seq import Seq


# read in the whole dataframe
df = pd.read_csv("../Data/Complete_ASOtoTranscriptSeq.tsv", sep=" ")

# list unique geneIDs
unique_genes = df.GeneID.unique()

# determine number of times each gene has aso in sequence
gene_check = {}
for i in unique_genes:
    gene_check[i] = 0

for i, r in df.iterrows():
    if str(Seq(r.ASOseq).reverse_complement()) in r.Sequence:
        gene_check[str(r.GeneID)] = gene_check[str(r.GeneID)] + 1

# identify genes with no matching ASOs to sequence
missing_genes = []
for i in gene_check.keys():
    if gene_check[i] == 0:
        missing_genes.append(i)


# TODO: need to figure out how we can have high ASOeffective but no ASO match in any of the transcripts?

# TODO: which genes do the missing ASOs belong to?

# TODO: check ASOs in blast to look for mismatch in ASOs?