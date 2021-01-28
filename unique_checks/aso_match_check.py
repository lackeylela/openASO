### checks for asos that never match any of the provided sequences

# import necessary modules
import pandas as pd
from Bio.Seq import Seq


# read in the whole dataframe
df = pd.read_csv("../Data/Complete_ASOtoTranscriptSeq.tsv", sep=" ")

# list unique geneIDs
unique_asos = df.ASOseq.unique()

# determine number of times each gene has aso in sequence
aso_check = {}
for i in unique_asos:
    aso_check[i] = 0

for i, r in df.iterrows():
    if str(Seq(r.ASOseq).reverse_complement()) in r.Sequence:
        aso_check[str(r.ASOseq)] = aso_check[str(r.ASOseq)] + 1

# identify genes with no matching ASOs to sequence
missing_asos = []
for i in aso_check.keys():
    if aso_check[i] == 0:
        missing_asos.append(i)

